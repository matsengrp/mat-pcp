#!/usr/bin/env python3
"""Pin the PCP gzip policy and the things it can silently break (#134).

The four PCP families are written as ``.csv.gz``. That is a rename plus a
compression policy, and three of its consequences are the kind that pass a
pipeline run while being wrong:

1. **The cross-implementation sentinel.** ``validate_test.py`` asserts the
   DASM-ready aggregate md5 equals ``dasm2.preprocess_command.preprocess``'s
   output, which is plain CSV. Compressing changes the file's bytes, so the md5
   has to be taken on the *decompressed* stream -- and then the expected value is
   unchanged. Re-baselining it to a compressed-file md5 would leave a green test
   that no longer compares anything to dasm2.
2. **Reproducibility.** pandas' default gzip output embeds ``mtime=now``, so an
   unchanged re-run would produce a different md5 -- which matters because
   ``data-central`` manifests pin ``md5``/``xxh3`` per file, and a spurious change
   reads as a data change.
3. **Clade labels.** ``aggregate_clade_pcps.py`` derives the ``clade`` column
   from each input's filename. Its old ``replace("_pcps.csv", "")`` would have
   yielded ``"19A.gz"`` under the new extension: every row of the aggregate
   mislabelled, identically, with the row count unchanged and no cross-check able
   to see it.

Checks 2 and 3 are mutation-tested -- each demonstrates that the failure it
guards is real -- because a check that cannot fail is worse than no check.

Usage:
    python tests/test_pcp_compression.py
"""

import gzip
import hashlib
import re
import sys
import tempfile
import traceback
from importlib import util as importlib_util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_DIR = REPO_ROOT / "workflow"

# The four families that must be gzipped, as they appear in Snakefile paths.
PCP_FAMILY_STEMS = (
    "{clade}_pcps",
    "all_pcps",
    "{clade}_pcps_dasm_ready",
    "all_pcps_dasm_ready",
)

# Deliberately *not* gzipped: single-column ID files, kept shell-readable.
PLAIN_CSV_STEMS = ("{clade}_filtered_branches", "filtered_branch_ids")


def _load_module(name, path):
    spec = importlib_util.spec_from_file_location(name, path)
    module = importlib_util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _load_pcpio():
    return _load_module("pcpio", WORKFLOW_DIR / "modules" / "pcpio.py")


def _exec_aggregate_script():
    """Execute ``aggregate_clade_pcps.py`` with a stub ``snakemake`` and return its globals.

    The script reads ``snakemake.input``/``.output``/``.log`` at module level, so
    it cannot be imported. ``__name__`` is set to something other than
    ``"__main__"`` so its ``main()`` does not run. This exercises the *shipped*
    ``clade_from_pcp_path`` rather than a copy of it.
    """
    path = WORKFLOW_DIR / "scripts" / "aggregate_clade_pcps.py"

    class _Stub:
        def __init__(self, **kw):
            self.__dict__.update(kw)

    stub = _Stub(
        input=_Stub(pcps=[]),
        output=_Stub(all_pcps="unused.csv.gz"),
        log=[],
    )
    namespace = {
        "__name__": "aggregate_clade_pcps_probe",
        "__file__": str(path),
        "snakemake": stub,
    }
    exec(compile(path.read_text(), str(path), "exec"), namespace)
    return namespace


def _sample_frame():
    """A frame shaped like a PCP table: long, near-identical sequences."""
    import pandas as pd

    parents = ["ACGT" * 60] * 40
    children = [("ACGT" * 60)[:-4] + "TTTT"] * 40
    return pd.DataFrame(
        {
            "clade": ["21J"] * 40,
            "parent_name": [f"node_{i}" for i in range(40)],
            "child_name": [f"node_{i + 1}" for i in range(40)],
            "parent": parents,
            "child": children,
            "branch_length": [1] * 40,
        }
    )


def check_policy_shape():
    """mtime=0 must be set, and compresslevel must be left at pandas' default."""
    failures = []
    pcpio = _load_pcpio()
    policy = pcpio.PCP_CSV_COMPRESSION

    if policy.get("method") != "gzip":
        failures.append(f"method is {policy.get('method')!r}, expected 'gzip'")
    if policy.get("mtime") != 0:
        failures.append(
            f"mtime is {policy.get('mtime')!r}, expected 0 -- without it an "
            f"unchanged re-run produces different bytes"
        )
    if "compresslevel" in policy:
        failures.append(
            f"compresslevel is pinned to {policy['compresslevel']!r}; leave it "
            f"unset so pandas' default (9) applies. Level 9 measured both "
            f"smaller and faster than 6 on PCP data (~125x vs ~17x on ORF1a)"
        )
    if not failures:
        print("  policy is gzip + mtime=0, compresslevel left at pandas' default")
    return failures


def check_round_trip():
    """A written .csv.gz must read back equal, and be genuinely compressed."""
    failures = []
    import pandas as pd

    pcpio = _load_pcpio()
    df = _sample_frame()

    with tempfile.TemporaryDirectory() as tmp:
        gz_path = Path(tmp) / "all_pcps.csv.gz"
        plain_path = Path(tmp) / "all_pcps.csv"
        df.to_csv(gz_path, index=False, compression=pcpio.PCP_CSV_COMPRESSION)
        df.to_csv(plain_path, index=False)

        back = pd.read_csv(gz_path)
        if not back.equals(df):
            failures.append("round-trip through .csv.gz did not preserve the frame")

        gz_size = gz_path.stat().st_size
        plain_size = plain_path.stat().st_size
        if gz_size >= plain_size:
            failures.append(
                f"gzip did not shrink the file ({gz_size} >= {plain_size})"
            )
        else:
            print(
                f"  round-trip exact; {plain_size} -> {gz_size} bytes "
                f"({plain_size / gz_size:.1f}x on this fixture)"
            )
    return failures


def check_decompressed_md5_matches_plain_csv():
    """The decompressed stream must equal the plain CSV byte for byte.

    This is the property that keeps ``DASM_READY_ALL_SENTINEL_MD5`` meaningful
    across the rename: our gzipped output decompresses to exactly what the
    plain-CSV reference implementation writes.
    """
    failures = []
    pcpio = _load_pcpio()
    df = _sample_frame()

    with tempfile.TemporaryDirectory() as tmp:
        gz_path = Path(tmp) / "all_pcps_dasm_ready.csv.gz"
        plain_path = Path(tmp) / "reference.csv"
        df.to_csv(gz_path, index=False, compression=pcpio.PCP_CSV_COMPRESSION)
        df.to_csv(plain_path, index=False)

        with gzip.open(gz_path, "rb") as handle:
            decompressed_md5 = hashlib.md5(handle.read()).hexdigest()
        plain_md5 = hashlib.md5(plain_path.read_bytes()).hexdigest()
        compressed_md5 = hashlib.md5(gz_path.read_bytes()).hexdigest()

        if decompressed_md5 != plain_md5:
            failures.append(
                f"decompressed md5 {decompressed_md5} != plain-CSV md5 {plain_md5}; "
                f"the cross-implementation sentinel would no longer hold"
            )
        # Mutation half: prove hashing the *container* would have been wrong.
        if compressed_md5 == plain_md5:
            failures.append(
                "compressed-file md5 equals the plain-CSV md5, so this check "
                "cannot distinguish the two and proves nothing"
            )
        if not failures:
            print("  decompressed stream is byte-identical to the plain CSV")
    return failures


def check_reruns_are_byte_identical():
    """Two writes of the same frame to the same path must match exactly.

    Mutation-tested: the same two writes *without* ``mtime=0`` must differ, or
    this check would pass for the wrong reason (e.g. a filesystem with
    second-granularity timestamps and a fast test).
    """
    failures = []
    import pandas as pd

    pcpio = _load_pcpio()
    df = _sample_frame()

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "all_pcps.csv.gz"

        # gzip's mtime field has 1-second resolution, so two back-to-back writes
        # can share a timestamp by luck. Force the difference by writing an
        # explicit non-zero mtime rather than sleeping.
        df.to_csv(path, index=False, compression={"method": "gzip", "mtime": 1})
        drifting_a = hashlib.md5(path.read_bytes()).hexdigest()
        df.to_csv(path, index=False, compression={"method": "gzip", "mtime": 2})
        drifting_b = hashlib.md5(path.read_bytes()).hexdigest()
        if drifting_a == drifting_b:
            failures.append(
                "gzip output did not change when mtime changed, so this check "
                "cannot detect a missing mtime=0"
            )

        df.to_csv(path, index=False, compression=pcpio.PCP_CSV_COMPRESSION)
        first = hashlib.md5(path.read_bytes()).hexdigest()
        df.to_csv(path, index=False, compression=pcpio.PCP_CSV_COMPRESSION)
        second = hashlib.md5(path.read_bytes()).hexdigest()
        if first != second:
            failures.append(
                f"two identical writes produced different bytes ({first} vs "
                f"{second}); mtime=0 is not taking effect"
            )
        elif not failures:
            print("  re-writing the same frame to the same path is byte-identical")
    return failures


def check_clade_derivation():
    """Clade names must survive the .csv.gz suffix."""
    failures = []
    namespace = _exec_aggregate_script()
    derive = namespace.get("clade_from_pcp_path")
    if derive is None:
        return [
            "aggregate_clade_pcps.py has no module-level clade_from_pcp_path; "
            "this check would pass vacuously"
        ]

    cases = {
        "results/x_v9/pcps/S/19A_pcps.csv.gz": "19A",
        "results/x_v9/pcps/ORF1a/21J_pcps.csv.gz": "21J",
        "24H_pcps.csv.gz": "24H",
        # The pre-#134 form must still resolve, for reading older runs.
        "19A_pcps.csv": "19A",
        # Clade names containing dots must not be truncated.
        "recombinant.1_pcps.csv.gz": "recombinant.1",
    }
    for path, expected in cases.items():
        got = derive(path)
        if got != expected:
            failures.append(f"{path} -> {got!r}, expected {expected!r}")

    # Mutation half: the old implementation must fail these, or the check is
    # not actually pinning the bug.
    import os as _os

    old_style = _os.path.basename("19A_pcps.csv.gz").replace("_pcps.csv", "")
    if old_style == "19A":
        failures.append(
            "the pre-fix replace() form yields '19A' here, so this check does "
            "not pin the #134 mislabelling bug"
        )

    try:
        derive("something_else.txt")
    except ValueError:
        pass
    else:
        failures.append("an unrecognized filename did not raise")

    if not failures:
        print("  clade derivation survives .csv.gz (old replace() would give '19A.gz')")
    return failures


def check_snakefile_paths():
    """Every PCP-family path in the Snakefile must end in .csv.gz.

    Guards the case where a future rule adds a fifth PCP output and writes it as
    plain CSV -- the outputs would be fine, but the directory would silently
    drift back and ``data-central``'s generator would emit a wrong path.
    """
    failures = []
    text = (WORKFLOW_DIR / "Snakefile").read_text()

    for stem in PCP_FAMILY_STEMS:
        bad = re.findall(rf'"[^"\n]*{re.escape(stem)}\.csv"', text)
        if bad:
            failures.append(
                f"{stem}: still referenced as plain .csv in the Snakefile: {bad}"
            )
        if not re.search(rf'{re.escape(stem)}\.csv\.gz"', text):
            failures.append(
                f"{stem}: no .csv.gz reference found in the Snakefile at all "
                f"(renamed or removed? this check may be stale)"
            )

    for stem in PLAIN_CSV_STEMS:
        if re.search(rf'{re.escape(stem)}\.csv\.gz"', text):
            failures.append(
                f"{stem}: gzipped, but it is documented as deliberately plain "
                f"CSV -- update OUTPUT_SCHEMA.md and this test together"
            )

    if not failures:
        print("  all four PCP families are .csv.gz; filtered-branch files stay .csv")
    return failures


def check_extension_is_dasm2_loadable():
    """Names must end in exactly '.csv.gz'.

    ``dasm2.data.registry.load_dataframe`` matches on the extension and accepts
    only ``.csv``/``.csv.gz`` (plus pickle variants), raising otherwise. dasm2 is
    not installed here, so this replicates that check rather than importing it.
    """
    failures = []
    accepted = (".csv", ".csv.gz", ".pkl", ".pickle", ".pkl.gz", ".pickle.gz")
    text = (WORKFLOW_DIR / "Snakefile").read_text()

    checked = 0
    for stem in PCP_FAMILY_STEMS:
        for match in re.findall(rf'"([^"\n]*{re.escape(stem)}[^"\n]*)"', text):
            # Only data paths matter here. A PCP stem also appears inside log and
            # sentinel filenames (e.g. ``log/curate_all_pcps_for_dasm_{gene}.log``),
            # which dasm2 never loads.
            if ".csv" not in match and ".gz" not in match:
                continue
            checked += 1
            if not match.lower().endswith(accepted):
                failures.append(
                    f"{match!r} would be rejected by dasm2's loader "
                    f"(accepted: {accepted}) -- note the check is on the exact "
                    f"extension, so '.gz' alone or '.csv.zst' does not work"
                )
    if not checked:
        failures.append(
            "no PCP data paths were examined; the stem list or the Snakefile "
            "layout has drifted and this check is passing vacuously"
        )
    if not failures:
        print(
            f"  all {checked} PCP data paths would be accepted by dasm2's "
            f"extension check"
        )
    return failures


CHECKS = [
    ("compression policy shape", check_policy_shape),
    ("round-trip through .csv.gz", check_round_trip),
    ("decompressed stream matches plain CSV", check_decompressed_md5_matches_plain_csv),
    ("re-runs are byte-identical", check_reruns_are_byte_identical),
    ("clade derivation from filename", check_clade_derivation),
    ("Snakefile PCP paths", check_snakefile_paths),
    ("dasm2 loader would accept the names", check_extension_is_dasm2_loadable),
]


def run_all_checks():
    failures = []
    for label, check in CHECKS:
        print(f"Checking {label}...")
        try:
            failures.extend(check())
        except Exception:
            failures.append(f"{label} raised:\n{traceback.format_exc()}")
    return failures


def test_pcp_compression():
    """Pytest entry point.

    This repo runs tests as standalone scripts (``pixi run python tests/<f>.py``)
    and does not depend on pytest. This wrapper exists so that anyone who does
    reach for ``pytest tests/`` gets a real result rather than a green run that
    collected nothing.
    """
    failures = run_all_checks()
    assert not failures, "\n".join(failures)


def main():
    failures = run_all_checks()
    if failures:
        print(f"\nFAILED ({len(failures)} problem(s)):")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("\nPASSED: PCP outputs are gzipped reproducibly, with payload-level md5s.")


if __name__ == "__main__":
    main()
