#!/usr/bin/env python3
"""Pin the founder-codon-background filter in the count table (#130).

``build_count_table.py`` joins observed counts onto the grid of every possible
mutation per clade. It used to join on ``(clade, nt_mutation)`` alone, dropping
``codon_change`` -- so an observation of the same nucleotide substitution on a
codon that had *already* mutated was relabelled onto the founder codon and kept
as a second row for the same mutation. In v8 that produced 6,415 duplicated
``(clade, nt_mutation)`` groups over 13,146 rows, affecting every clade and
concentrating in the highest-count mutations (62.5% of 21J spike mutations with
>=1000 counts).

Both fitness estimators were corrupted, in opposite directions and by different
mechanisms, so neither a schema check nor a cross-check between them could see
it -- and the amino-acid uncertainty got *smaller*, so it could not be filtered
out downstream either.

What this pins:

* **the filter** -- an off-founder-background observation is dropped, and the
  surviving row carries the founder-background count only;
* **the paths that must survive it** -- noncoding mutations (no codon at all)
  and overlap sites whose ``;``-joined codons differ between genes, which is
  where a naive comparison would silently zero a real count;
* **the row-count invariant** -- one row per ``(clade, nt_mutation)``, which is
  the guarantee the fused rule lost;
* **the structural guard** -- a coding/gene-structure disagreement between
  ``codon_change`` and the founder grid raises rather than being read as
  "off-founder";
* **the drop log** -- what was removed is reported per clade, since mutations
  recurring on mutated backgrounds are a real quantity, not noise;
* **the independent guard** in ``aamutfit.aa_fitness``, which would have caught
  this on its own.

The test-data pipeline cannot cover any of this: its 66 observed rows contain
zero duplicated keys and zero off-founder backgrounds, so a green test-data run
says nothing about the filter. Hence synthetic fixtures.

``build_count_table.py`` is a Snakemake script driven by a module-level
``snakemake`` global, so -- as in ``test_clade_founder_guards.py`` -- it is run
the way Snakemake runs it, with a stub injected. That keeps the shipped code
path under test rather than a testable extract of it.

Usage:
    python tests/test_founder_codon_background.py
"""

import csv
import sys
import tempfile
import traceback
import types
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "build_count_table.py"


def _ensure_biopython():
    """Make ``import Bio.SeqIO`` work, with a stub if biopython is absent.

    biopython lives in ``workflow/envs/default.yml``, not the pixi environment:
    this repo deliberately keeps science dependencies in the conda envs. Rather
    than widen the pixi env for one test, fall back to a minimal single-record
    FASTA reader. Real biopython is preferred when present, so a conda-env run
    exercises the actual parser.
    """
    try:
        import Bio.SeqIO  # noqa: F401

        return "biopython"
    except ModuleNotFoundError:
        pass

    class _Record:
        def __init__(self, seq):
            self.seq = seq

    def _read(path, fmt):
        assert fmt == "fasta", f"stub only handles fasta, got {fmt!r}"
        lines = Path(path).read_text().splitlines()
        # Reject a non-header first line rather than skipping it: `lines[1:]`
        # assumes line 0 is the header, so without this a file that opens with
        # sequence data silently loses its first line. This check exists in
        # tests/test_clade_founder_guards.py's copy of this stub and had drifted
        # out of this one -- see the PR #132 review.
        if not lines or not lines[0].startswith(">"):
            raise ValueError(f"{path} is not a FASTA file")
        body = [ln.strip() for ln in lines[1:] if ln.strip() and not ln.startswith(">")]
        return _Record("".join(body))

    bio = types.ModuleType("Bio")
    seqio = types.ModuleType("Bio.SeqIO")
    seqio.read = _read
    bio.SeqIO = seqio
    sys.modules.setdefault("Bio", bio)
    sys.modules.setdefault("Bio.SeqIO", seqio)
    return "stub"


# A founder with five sites, chosen to exercise every path through the filter:
#
#   site 10  codon GCT at S:67 pos 2   -- the A67V case from #130
#   site 20  noncoding                 -- no codon to compare
#   site 30  N;ORF9b overlap, codons AAA;AAC in different frames
#   site 40  codon TCG at S:221 pos 2  -- the 4-row S221L case
#   site 50  ORF1a;ORF1ab overlap, same frame so the codons are identical
#
# Sites are spaced so the reference sequence below can be short.
FOUNDER_NTS = [
    # clade, site, nt, gene, codon, codon_position, codon_site, four_fold_degenerate
    ("21J", 10, "C", "S", "GCT", 2, 67, False),
    ("21J", 20, "A", "noncoding", "noncoding", "noncoding", "noncoding", False),
    ("21J", 30, "A", "N;ORF9b", "AAA;AAC", "2;1", "38;35", False),
    ("21J", 40, "C", "S", "TCG", 2, 221, False),
    ("21J", 50, "C", "ORF1a;ORF1ab", "GCT;GCT", "2;2", "76;76", False),
]

REF_SEQ = "ACGT" * 15  # 60 nt; only the five sites above are ever looked up.


def _write_founder_nts(path, rows=None):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            ["clade", "site", "nt", "gene", "codon", "codon_position",
             "codon_site", "four_fold_degenerate"]
        )
        w.writerows(FOUNDER_NTS if rows is None else rows)


ACTUAL_HEADER = [
    "clade", "protein", "aa_mutation", "nt_mutation", "codon_change",
    "synonymous", "noncoding", "nt_site", "reference_nt", "clade_founder_nt",
    "count", "count_terminal", "count_non_terminal", "mean_log_size",
]


def _actual_row(clade, nt_mutation, codon_change, count, protein="S",
                aa_mutation="A67V", noncoding=False):
    return [
        clade, protein, aa_mutation, nt_mutation, codon_change, False, noncoding,
        int(nt_mutation[1:-1]), nt_mutation[0], nt_mutation[0], count, count, 0, 0.0,
    ]


def _write_actual(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(ACTUAL_HEADER)
        w.writerows(rows)


class _Namespace:
    """Stand-in for Snakemake's input/output/params objects."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def get(self, key, default=None):
        return self.__dict__.get(key, default)


def _run_script(tmpdir, actual_rows, founder_rows=None):
    """Execute build_count_table.py the way Snakemake does.

    Returns ``(rows, log_text)`` where rows is a list of dicts from the written
    ``counts_all.csv``. Raises whatever the script raises.
    """
    tmpdir = Path(tempfile.mkdtemp(dir=tmpdir))
    founder_csv = tmpdir / "clade_founder_nts.csv"
    actual_csv = tmpdir / "aggregated.csv"
    ref_fasta = tmpdir / "ref.fa"
    site_mask_csv = tmpdir / "site_mask.csv"
    usher_yaml = tmpdir / "usher_masked_sites.yaml"
    out_csv = tmpdir / "counts_all.csv"
    log_path = tmpdir / "build.log"

    _write_founder_nts(founder_csv, founder_rows)
    _write_actual(actual_csv, actual_rows)
    ref_fasta.write_text(f">ref\n{REF_SEQ}\n")
    site_mask_csv.write_text("site\n")
    usher_yaml.write_text("{}\n")

    stub = _Namespace(
        input=_Namespace(
            actual=str(actual_csv),
            clade_founder_nts_csv=str(founder_csv),
            ref_fasta=str(ref_fasta),
            site_mask=str(site_mask_csv),
            usher_masked_sites=str(usher_yaml),
        ),
        output=_Namespace(csv=str(out_csv)),
        params=_Namespace(
            sites_to_exclude=[],
            # Widened to the fixture genome so nothing is flagged
            # excluded; the filter under test runs before masking anyway.
            site_include_range={"start": 1, "end": len(REF_SEQ)},
        ),
        log=[str(log_path)],
    )

    script_globals = {
        "__name__": "__main__",
        "__file__": str(SCRIPT),
        "snakemake": stub,
    }
    exec(compile(SCRIPT.read_text(), str(SCRIPT), "exec"), script_globals)

    with open(out_csv) as f:
        rows = list(csv.DictReader(f))
    return rows, log_path.read_text()


def _find(rows, clade, nt_mutation):
    return [r for r in rows if r["clade"] == clade and r["nt_mutation"] == nt_mutation]


def check_off_founder_background_dropped(tmpdir):
    """The #130 worked example: A67V keeps 2254, not 2254 and a phantom 1."""
    failures = []
    rows, log_text = _run_script(
        tmpdir,
        [
            _actual_row("21J", "C10T", "GCT>GTT", 2254, aa_mutation="A67V"),
            # Same nucleotide event on a codon that had already taken A67S.
            # This is S67F, not A67V, and must not contribute to 21J A67V.
            _actual_row("21J", "C10T", "TCT>TTT", 1, aa_mutation="S67F"),
        ],
    )

    matched = _find(rows, "21J", "C10T")
    if len(matched) != 1:
        return [
            f"off-founder background: expected exactly 1 row for (21J, C10T), "
            f"got {len(matched)}. This is the #130 failure mode -- the phantom "
            "row would get ~1/3 of the inverse-variance weight."
        ]
    if int(matched[0]["actual_count"]) != 2254:
        failures.append(
            f"off-founder background: actual_count is "
            f"{matched[0]['actual_count']}, expected 2254 (the "
            "founder-background count alone, neither summed with nor replaced "
            "by the off-founder observation)"
        )
    if matched[0]["aa_mutation"] != "A67V":
        failures.append(
            f"off-founder background: aa_mutation is "
            f"{matched[0]['aa_mutation']!r}, expected 'A67V'"
        )
    # The dropped observation is a real biological quantity; it must be
    # reported, not silently discarded.
    if "mutated codon background" not in log_text or "21J" not in log_text:
        failures.append(
            "off-founder background: the drop was not logged per clade. "
            f"Log was:\n{log_text}"
        )
    if not failures:
        print("  off-founder codon background: dropped, count preserved, logged")
    return failures


def check_multiple_off_founder_rows(tmpdir):
    """The S221L case: 4 observed rows collapse to 1, keeping only the founder's."""
    failures = []
    rows, _ = _run_script(
        tmpdir,
        [
            _actual_row("21J", "C40T", "TCG>TTG", 1218, aa_mutation="S221L"),
            _actual_row("21J", "C40T", "GCG>GTG", 1, aa_mutation="A221V"),
            _actual_row("21J", "C40T", "TCT>TTT", 1, aa_mutation="S221F"),
            _actual_row("21J", "C40T", "TCA>TTA", 1, aa_mutation="S221L"),
        ],
    )
    matched = _find(rows, "21J", "C40T")
    if len(matched) != 1:
        failures.append(
            f"multi-background: expected 1 row for (21J, C40T), got "
            f"{len(matched)}. Upstream reported expected_count 4x too large "
            "here, turning a near-neutral mutation into a strongly deleterious "
            "one."
        )
    elif int(matched[0]["actual_count"]) != 1218:
        failures.append(
            f"multi-background: actual_count is {matched[0]['actual_count']}, "
            "expected 1218"
        )
    if not failures:
        print("  four codon backgrounds: collapsed to the founder row only")
    return failures


def check_noncoding_and_overlap_survive(tmpdir):
    """Noncoding and out-of-frame overlap observations must pass the filter.

    These are the two ways a correct observation could be mistaken for an
    off-founder one and silently zeroed -- which the row-count invariant would
    *not* catch, since the grid row survives either way.
    """
    failures = []
    rows, _ = _run_script(
        tmpdir,
        [
            _actual_row("21J", "A20G", "noncoding", 7, protein="noncoding",
                        aa_mutation="noncoding", noncoding=True),
            # Overlap site: the two genes are in different frames, so the
            # ;-joined codons genuinely differ. Comparing them requires the
            # same gene order on both sides.
            _actual_row("21J", "A30T", "AAA>ATA;AAC>TAC", 5,
                        protein="N;ORF9b", aa_mutation="K38I;N35Y"),
            # ORF1a/ORF1ab: the same frame, so the ;-joined codons are
            # identical. This is the weaker of the two overlap cases -- it
            # would pass even with a broken comparison -- but it is the
            # overlap that actually carries most of the genome, so the
            # ;-join path is pinned for it too.
            _actual_row("21J", "C50T", "GCT>GTT;GCT>GTT", 9,
                        protein="ORF1a;ORF1ab", aa_mutation="A76V;A76V"),
        ],
    )
    for nt_mutation, expected, label in (
        ("A20G", 7, "noncoding"),
        ("A30T", 5, "out-of-frame overlap"),
        ("C50T", 9, "same-frame overlap"),
    ):
        matched = _find(rows, "21J", nt_mutation)
        if len(matched) != 1:
            failures.append(
                f"{label}: expected 1 row for {nt_mutation}, got {len(matched)}"
            )
        elif int(matched[0]["actual_count"]) != expected:
            failures.append(
                f"{label}: {nt_mutation} kept actual_count "
                f"{matched[0]['actual_count']}, expected {expected} -- a "
                "legitimate observation was dropped to zero by the filter"
            )
    if not failures:
        print("  noncoding and out-of-frame overlap observations survive")
    return failures


def check_unobserved_rows_stay_zero(tmpdir):
    """Every grid row survives, unobserved ones at zero."""
    failures = []
    rows, _ = _run_script(
        tmpdir, [_actual_row("21J", "C10T", "GCT>GTT", 3, aa_mutation="A67V")]
    )
    expected_rows = len(FOUNDER_NTS) * 3  # 3 alternative bases per site
    if len(rows) != expected_rows:
        failures.append(
            f"grid invariant: wrote {len(rows)} rows, expected {expected_rows} "
            f"({len(FOUNDER_NTS)} sites x 3 alternative bases). The grid is the "
            "authoritative enumeration and must survive the join unchanged."
        )
    observed = [r for r in rows if int(r["actual_count"]) > 0]
    if len(observed) != 1:
        failures.append(
            f"grid invariant: {len(observed)} rows have a nonzero count, "
            "expected 1"
        )
    if not failures:
        print("  full grid preserved, unobserved mutations at zero")
    return failures


def check_structural_mismatch_raises(tmpdir):
    """A coding/gene-structure disagreement raises instead of dropping."""
    try:
        _run_script(
            tmpdir,
            [
                # The grid says site 30 is covered by two genes in different
                # frames; this observation names only one codon. Reading that
                # as "off-founder" would silently zero the count.
                _actual_row("21J", "A30T", "AAA>ATA", 5, protein="N",
                            aa_mutation="K38I"),
            ],
        )
    except ValueError as exc:
        message = str(exc)
        missing = [
            m for m in ("coding/gene structure", "#130", "A30T") if m not in message
        ]
        if missing:
            return [
                f"structural mismatch: raised, but the message omits {missing}. "
                f"Got:\n{message}"
            ]
        print("  structural mismatch: raised ValueError as required")
        return []
    except Exception:  # noqa: BLE001
        return [f"structural mismatch: wrong exception:\n{traceback.format_exc()}"]
    return [
        "structural mismatch: did NOT raise. A gene-structure disagreement "
        "would be read as an off-founder observation and silently zero a real "
        "count -- and the row-count invariant cannot see it."
    ]


def check_missing_codon_change_raises(tmpdir):
    """A counts file without codon_change is a hard error, not a skipped filter."""
    tmp = Path(tempfile.mkdtemp(dir=tmpdir))
    actual_csv = tmp / "aggregated.csv"
    # Same rows, codon_change column removed entirely.
    keep = [i for i, c in enumerate(ACTUAL_HEADER) if c != "codon_change"]
    row = _actual_row("21J", "C10T", "GCT>GTT", 3)
    with open(actual_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([ACTUAL_HEADER[i] for i in keep])
        w.writerow([row[i] for i in keep])

    founder_csv = tmp / "clade_founder_nts.csv"
    _write_founder_nts(founder_csv)
    ref_fasta = tmp / "ref.fa"
    ref_fasta.write_text(f">ref\n{REF_SEQ}\n")
    (tmp / "site_mask.csv").write_text("site\n")
    (tmp / "usher.yaml").write_text("{}\n")

    stub = _Namespace(
        input=_Namespace(
            actual=str(actual_csv),
            clade_founder_nts_csv=str(founder_csv),
            ref_fasta=str(ref_fasta),
            site_mask=str(tmp / "site_mask.csv"),
            usher_masked_sites=str(tmp / "usher.yaml"),
        ),
        output=_Namespace(csv=str(tmp / "counts_all.csv")),
        params=_Namespace(
            sites_to_exclude=[],
            # Widened to the fixture genome so nothing is flagged
            # excluded; the filter under test runs before masking anyway.
            site_include_range={"start": 1, "end": len(REF_SEQ)},
        ),
        log=[str(tmp / "build.log")],
    )
    try:
        exec(
            compile(SCRIPT.read_text(), str(SCRIPT), "exec"),
            {"__name__": "__main__", "__file__": str(SCRIPT), "snakemake": stub},
        )
    except ValueError as exc:
        if "codon_change" not in str(exc):
            return [f"missing codon_change: raised the wrong ValueError: {exc}"]
        print("  missing codon_change column: raised ValueError as required")
        return []
    except Exception:  # noqa: BLE001
        return [f"missing codon_change: wrong exception:\n{traceback.format_exc()}"]
    return [
        "missing codon_change: did NOT raise. Without that column the filter "
        "cannot run at all, and skipping it silently is exactly how #130 stayed "
        "invisible."
    ]


def check_aa_fitness_duplicate_guard():
    """``aamutfit.aa_fitness`` rejects a repeated nt_mutation independently."""
    sys.path.insert(0, str(REPO_ROOT / "workflow" / "modules"))
    try:
        import pandas as pd

        import aamutfit
    except ModuleNotFoundError as exc:
        return [f"aa_fitness guard: could not import ({exc})"]

    # Two rows for the same nucleotide mutation in one amino-acid group: the
    # shape #130 produced downstream. Averaging them gives a phantom
    # observation real weight and a *smaller* uncertainty.
    df = pd.DataFrame(
        {
            "cluster": ["21J", "21J"],
            "gene": ["S", "S"],
            "codon_site": [67, 67],
            "aa_mutation": ["A67V", "A67V"],
            "nt_mutation": ["C21762T", "C21762T"],
            "f_mean": [2.823, -4.101],
            "f_st_dev": [0.695, 1.036],
            "actual_count": [2254, 1],
        }
    )
    explode_cols = ["gene", "codon_site", "aa_mutation"]
    try:
        aamutfit.aa_fitness(df, explode_cols)
    except ValueError as exc:
        if "#130" not in str(exc) or "nt_mutation" not in str(exc):
            return [f"aa_fitness guard: raised, but message is unhelpful: {exc}"]
        print("  aa_fitness duplicate guard: raised ValueError as required")
    except Exception:  # noqa: BLE001
        return [f"aa_fitness guard: wrong exception:\n{traceback.format_exc()}"]
    else:
        return [
            "aa_fitness guard: did NOT raise on a repeated nt_mutation within "
            "one amino-acid group. This is the cheapest independent check on "
            "#130 and would have caught it on its own."
        ]

    # And it must not fire on genuinely distinct nucleotide paths to the same
    # amino-acid change, which is the case the average exists for.
    df.loc[1, "nt_mutation"] = "C21763T"
    try:
        out = aamutfit.aa_fitness(df, explode_cols)
    except Exception:  # noqa: BLE001
        return [
            "aa_fitness guard: fired on two *distinct* nucleotide paths to the "
            f"same amino-acid change, which is legitimate:\n{traceback.format_exc()}"
        ]
    if len(out) != 1:
        return [f"aa_fitness guard: expected 1 averaged row, got {len(out)}"]
    print("  aa_fitness: distinct nucleotide paths still averaged")
    return []


def run_all_checks():
    failures = []
    if not SCRIPT.is_file():
        return [f"missing script {SCRIPT}"]

    print(f"  FASTA reader: {_ensure_biopython()}")

    with tempfile.TemporaryDirectory() as tmpdir:
        print("Checking off-founder codon backgrounds are dropped (#130)...")
        failures += check_off_founder_background_dropped(tmpdir)
        failures += check_multiple_off_founder_rows(tmpdir)

        print("Checking legitimate observations survive the filter...")
        failures += check_noncoding_and_overlap_survive(tmpdir)
        failures += check_unobserved_rows_stay_zero(tmpdir)

        print("Checking the guards fail loudly...")
        failures += check_structural_mismatch_raises(tmpdir)
        failures += check_missing_codon_change_raises(tmpdir)

    print("Checking the independent aa-level guard...")
    failures += check_aa_fitness_duplicate_guard()

    return failures


def test_founder_codon_background():
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
    print(
        "\nPASSED: only founder-codon-background observations are counted, "
        "one row per mutation."
    )


if __name__ == "__main__":
    main()
