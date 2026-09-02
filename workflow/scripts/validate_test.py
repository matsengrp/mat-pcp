#!/usr/bin/env python3
"""Validate test pipeline outputs against expected results."""

import json
import math
import os
import re
import sys
from pathlib import Path

# Cross-version integrity sentinel for the gene-sliced 19A spike founder.
# 19A spike is byte-stable across SARS-CoV-2 mat-pcp cutoff dates because the
# Wuhan-Hu-1 reference founder for clade 19A doesn't change between runs.
# A mismatch here means either the slice math drifted, the reference FASTA
# changed, or the FASTA writer's header/wrap convention changed.
SPIKE_19A_SENTINEL_MD5 = "b948f2a1f20edbce64cfff1dc540c119"

# Cross-implementation byte-equivalence sentinels for DASM-ready PCP outputs
# on the test dataset (#80). These md5s come from running the reference
# implementation in ``dasm2.preprocess_command.preprocess`` (the same
# function that ``data-central/scripts/filter_viral_pcps.sh`` invokes) on
# the test pipeline's ``all_pcps.csv`` + ``S_codon_sites_to_mask.txt``
# with the canonical seeds. A mismatch here means either the filter logic
# drifted from dasm2 or the train/val RNG call sequence changed.
#
# Since #134 our side is gzipped while the reference writes plain CSV, so
# ``test_summary.py`` takes these md5s over the **decompressed** stream
# (``_decompressed_md5``). The sentinel value is therefore *unchanged* by
# compression -- which is the point: it still compares payload to payload. Do
# not re-baseline it to a compressed-file md5. That would silently retire the
# cross-implementation check, and it would not even be stable, since pandas
# embeds the output filename in the gzip header (FNAME) and, unless
# ``mtime=0`` is passed, the write time as well.
DASM_READY_ALL_SENTINEL_MD5 = "5e538bd088e676b485a11a04745221e9"


def is_close(actual, expected, rel_tol=1e-6, abs_tol=1e-9):
    """Check if two floats are close within tolerance."""
    if expected == 0:
        return abs(actual) <= abs_tol
    return abs(actual - expected) <= max(rel_tol * abs(expected), abs_tol)


def load_json(filepath):
    """Load JSON from file."""
    with open(filepath, 'r') as f:
        return json.load(f)


def compare_summaries(actual, expected):
    """Compare actual and expected summaries, returning differences."""
    differences = []

    # Check total PCPs
    if actual['pcp_summary']['total_pcps'] != expected['pcp_summary']['total_pcps']:
        differences.append(
            f"PCP count mismatch: expected {expected['pcp_summary']['total_pcps']}, "
            f"got {actual['pcp_summary']['total_pcps']}"
        )

    # Check columns
    actual_cols = set(actual['pcp_summary']['columns'])
    expected_cols = set(expected['pcp_summary']['columns'])
    if actual_cols != expected_cols:
        differences.append(
            f"Column mismatch: missing {expected_cols - actual_cols}, "
            f"extra {actual_cols - expected_cols}"
        )

    # Check clade counts
    if actual['clade_counts_summary']['total_clades'] != expected['clade_counts_summary']['total_clades']:
        differences.append(
            f"Clade count mismatch: expected {expected['clade_counts_summary']['total_clades']}, "
            f"got {actual['clade_counts_summary']['total_clades']}"
        )

    # Check total sequences
    if actual['clade_counts_summary']['total_sequences'] != expected['clade_counts_summary']['total_sequences']:
        differences.append(
            f"Total sequence count mismatch: expected {expected['clade_counts_summary']['total_sequences']}, "
            f"got {actual['clade_counts_summary']['total_sequences']}"
        )

    # Check clade counts checksum
    if actual['clade_counts_summary']['checksum'] != expected['clade_counts_summary']['checksum']:
        differences.append("Clade counts file has changed (checksum mismatch)")

    # Check for empty clade files
    if len(actual.get('empty_clade_files', [])) > 0:
        differences.append(f"Found empty clade files: {actual['empty_clade_files']}")

    # Check file existence (normalize paths by removing directory prefix)
    def normalize_path(path):
        """Strip the results directory prefix to compare relative paths."""
        # Remove results/test_data or results/test_data_v1 etc prefix
        return re.sub(r'^results/test_data[^/]*/', '', path)

    expected_normalized = {normalize_path(k): v for k, v in expected['files'].items()}
    actual_normalized = {normalize_path(k): v for k, v in actual['files'].items()}

    for filepath_norm, stats in expected_normalized.items():
        if filepath_norm not in actual_normalized:
            differences.append(f"Missing expected file: {filepath_norm}")
        elif actual_normalized[filepath_norm] is None:
            differences.append(f"Expected file does not exist: {filepath_norm}")

    # Check rates_table.csv
    if expected.get('rates_table_summary'):
        if not actual.get('rates_table_summary'):
            differences.append("rates_table.csv summary missing")
        else:
            exp_rates = expected['rates_table_summary']
            act_rates = actual['rates_table_summary']

            if act_rates['num_rows'] != exp_rates['num_rows']:
                differences.append(
                    f"rates_table row count mismatch: expected {exp_rates['num_rows']}, "
                    f"got {act_rates['num_rows']}"
                )
            if set(act_rates.get('mut_types', [])) != set(exp_rates.get('mut_types', [])):
                differences.append(
                    f"rates_table mutation types mismatch: expected {exp_rates['mut_types']}, "
                    f"got {act_rates['mut_types']}"
                )
            if act_rates.get('num_motifs') != exp_rates.get('num_motifs'):
                differences.append(
                    f"rates_table motif count mismatch: expected {exp_rates['num_motifs']}, "
                    f"got {act_rates['num_motifs']}"
                )
            # Check rate_sum with floating-point tolerance
            if 'rate_sum' in exp_rates and 'rate_sum' in act_rates:
                if not is_close(act_rates['rate_sum'], exp_rates['rate_sum'], rel_tol=1e-6):
                    differences.append(
                        f"rates_table rate_sum mismatch: expected {exp_rates['rate_sum']}, "
                        f"got {act_rates['rate_sum']}"
                    )

    # Check clade_founder_nts.csv
    if expected.get('clade_founder_nts_summary'):
        if not actual.get('clade_founder_nts_summary'):
            differences.append("clade_founder_nts.csv summary missing")
        else:
            exp_nts = expected['clade_founder_nts_summary']
            act_nts = actual['clade_founder_nts_summary']

            if act_nts['num_rows'] != exp_nts['num_rows']:
                differences.append(
                    f"clade_founder_nts row count mismatch: expected {exp_nts['num_rows']}, "
                    f"got {act_nts['num_rows']}"
                )
            if set(act_nts.get('clades', [])) != set(exp_nts.get('clades', [])):
                differences.append(
                    f"clade_founder_nts clades mismatch: expected {exp_nts['clades']}, "
                    f"got {act_nts['clades']}"
                )
            if act_nts['checksum'] != exp_nts['checksum']:
                differences.append("clade_founder_nts.csv has changed (checksum mismatch)")

    # Check gene-sliced (spike) founder FASTAs, including the cross-version
    # 19A integrity sentinel (see module-level ``SPIKE_19A_SENTINEL_MD5``).
    if expected.get('spike_founders_summary'):
        if not actual.get('spike_founders_summary'):
            differences.append("spike_founders_summary missing")
        else:
            exp_spike = expected['spike_founders_summary']
            act_spike = actual['spike_founders_summary']

            if act_spike.get('num_files') != exp_spike.get('num_files'):
                differences.append(
                    f"spike_founders file count mismatch: expected "
                    f"{exp_spike.get('num_files')}, got {act_spike.get('num_files')}"
                )
            if set(act_spike.get('files', [])) != set(exp_spike.get('files', [])):
                differences.append(
                    f"spike_founders files mismatch: expected "
                    f"{exp_spike.get('files')}, got {act_spike.get('files')}"
                )
            for clade, exp_md5 in exp_spike.get('checksums', {}).items():
                act_md5 = act_spike.get('checksums', {}).get(clade)
                if act_md5 != exp_md5:
                    differences.append(
                        f"spike_founders {clade} md5 mismatch: expected "
                        f"{exp_md5}, got {act_md5}"
                    )

            # Hard-coded cross-version sentinel: 19A spike slice must match
            # the canonical Wuhan-Hu-1 reference checksum.
            act_19a = act_spike.get('checksums', {}).get('19A')
            if act_19a and act_19a != SPIKE_19A_SENTINEL_MD5:
                differences.append(
                    f"19A spike sentinel md5 mismatch: expected "
                    f"{SPIKE_19A_SENTINEL_MD5}, got {act_19a}"
                )

    # Check DASM-ready PCP outputs (#80), including byte-equivalence sentinel
    # against the dasm2 reference implementation.
    if expected.get('dasm_ready_pcps_summary'):
        if not actual.get('dasm_ready_pcps_summary'):
            differences.append("dasm_ready_pcps_summary missing")
        else:
            exp_dasm = expected['dasm_ready_pcps_summary']
            act_dasm = actual['dasm_ready_pcps_summary']

            exp_agg = exp_dasm.get('aggregate', {})
            act_agg = act_dasm.get('aggregate', {})

            if act_agg.get('num_rows') != exp_agg.get('num_rows'):
                differences.append(
                    f"dasm_ready aggregate row count mismatch: expected "
                    f"{exp_agg.get('num_rows')}, got {act_agg.get('num_rows')}"
                )
            if not act_agg.get('has_in_train'):
                differences.append("dasm_ready aggregate missing in_train column")
            if act_agg.get('n_train') != exp_agg.get('n_train'):
                differences.append(
                    f"dasm_ready aggregate n_train mismatch: expected "
                    f"{exp_agg.get('n_train')}, got {act_agg.get('n_train')}"
                )
            # n_val asserted separately so a split bug that flips
            # membership without changing total count is caught even if
            # row count and n_train both happen to match.
            if act_agg.get('n_val') != exp_agg.get('n_val'):
                differences.append(
                    f"dasm_ready aggregate n_val mismatch: expected "
                    f"{exp_agg.get('n_val')}, got {act_agg.get('n_val')}"
                )
            if act_agg.get('checksum') != exp_agg.get('checksum'):
                differences.append(
                    f"dasm_ready aggregate md5 mismatch: expected "
                    f"{exp_agg.get('checksum')}, got {act_agg.get('checksum')}"
                )

            # Hard-coded cross-implementation sentinel: the aggregate md5
            # must equal the dasm2 reference output (proves byte-equivalence
            # to filter_viral_pcps.sh). Taken on the decompressed stream since
            # #134, so the expected value is unchanged by compression.
            if act_agg.get('checksum') and act_agg.get('checksum') != DASM_READY_ALL_SENTINEL_MD5:
                differences.append(
                    f"dasm_ready aggregate sentinel md5 mismatch: expected "
                    f"{DASM_READY_ALL_SENTINEL_MD5}, got {act_agg.get('checksum')}"
                )

            # Per-clade checksum diffs
            for clade, exp_md5 in exp_dasm.get('per_clade_checksums', {}).items():
                act_md5 = act_dasm.get('per_clade_checksums', {}).get(clade)
                if act_md5 != exp_md5:
                    differences.append(
                        f"dasm_ready per-clade {clade} md5 mismatch: expected "
                        f"{exp_md5}, got {act_md5}"
                    )

    # Check aamut_fitness outputs
    if expected.get('aamut_fitness_summary'):
        if not actual.get('aamut_fitness_summary'):
            differences.append("aamut_fitness summary missing")
        else:
            exp_aa = expected['aamut_fitness_summary']
            act_aa = actual['aamut_fitness_summary']

            if act_aa['num_rows'] != exp_aa['num_rows']:
                differences.append(
                    f"aamut_fitness row count mismatch: expected {exp_aa['num_rows']}, "
                    f"got {act_aa['num_rows']}"
                )
            if act_aa['num_clade_files'] != exp_aa['num_clade_files']:
                differences.append(
                    f"aamut_fitness clade file count mismatch: expected {exp_aa['num_clade_files']}, "
                    f"got {act_aa['num_clade_files']}"
                )
            if set(act_aa.get('clade_files', [])) != set(exp_aa.get('clade_files', [])):
                differences.append(
                    f"aamut_fitness clade files mismatch: expected {exp_aa['clade_files']}, "
                    f"got {act_aa['clade_files']}"
                )

            # Check numerical sums with appropriate tolerances
            # Integer-like sums (actual_count) use tight tolerance
            for key in ['actual_count_sum']:
                if key in exp_aa and key in act_aa:
                    if not is_close(act_aa[key], exp_aa[key], rel_tol=1e-6):
                        differences.append(
                            f"aamut_fitness {key} mismatch: expected {exp_aa[key]}, "
                            f"got {act_aa[key]}"
                        )
            # Floating-point sums can vary more across platforms, use looser tolerance
            for key in ['expected_count_sum', 'delta_fitness_ratio_sum',
                        'delta_fitness_bayesian_sum']:
                if key in exp_aa and key in act_aa:
                    if not is_close(act_aa[key], exp_aa[key], rel_tol=0.01):
                        rel_diff = (
                            abs(act_aa[key] - exp_aa[key]) / abs(exp_aa[key])
                            if exp_aa[key] else float('inf')
                        )
                        differences.append(
                            f"aamut_fitness {key} mismatch: expected {exp_aa[key]}, "
                            f"got {act_aa[key]} (rel_diff={rel_diff:.4f})"
                        )

    return differences


def validate_root_to_tip_cloud(output_dir):
    """The RTT cloud, the report row and the figure all describe the same regression.

    The scatter (#146) exists because ``chronumental_qc``'s aggregate check divides
    by chronumental's pre-fit root-to-tip slope, and that slope is the only
    *independent* reference in the QC. Three things a run can check that the unit
    tests cannot, because they are cross-*file* properties:

    1. ``node_dates_qc.tsv``'s ``root_to_tip_regression`` slope equals a
       recomputation from ``node_dates_rtt.csv``. The figure recomputes from that
       CSV too, so this is what makes "the figure cannot disagree with the report"
       an assertion rather than a design intention.
    2. The cloud's tips reconcile with the tips in ``node_dates_residuals.csv``.
       Both are the day-resolution set, and the cloud additionally requires the
       strain to be a node in *this* tree -- so the cloud may be smaller, never
       larger. A silently *larger* cloud would mean the strain->node mapping is
       matching something it should not.
    3. Divergence is non-negative and its own maximum is plausible. A negative
       value means the accumulation ran the wrong way down the tree, which the
       scatter would render as a mirrored cloud that still looks like a cloud.

    Deliberately **not** asserted: that the recomputed slope equals chronumental's
    logged one. The two regressions run over different tip sets (see
    ``node_date_qc.root_to_tip_histogram``), so equality is false in general -- the
    report states both and flags only a gross gap. It happens to hold to 0.1% on the
    test data, where every tip is day-resolution, and reading that as a guarantee is
    what would make this check fail on the production trees.
    """
    import pandas as pd

    cloud_path = os.path.join(output_dir, "mat", "node_dates_rtt.csv")
    report_path = os.path.join(output_dir, "mat", "node_dates_qc.tsv")
    plot_path = os.path.join(output_dir, "mat", "node_dates_rtt_scatter.png")
    if not os.path.exists(cloud_path):
        return [], f"SKIPPED: no {cloud_path}"

    differences = []
    cloud = pd.read_csv(cloud_path, parse_dates=["observed_date"])
    if not os.path.exists(plot_path) or os.path.getsize(plot_path) == 0:
        differences.append(
            f"root-to-tip: {plot_path} is missing or empty, so the cloud was "
            f"computed but never rendered"
        )

    if cloud.empty:
        # Legitimate only if there were no day-resolution tips at all.
        return differences, "root-to-tip: cloud is empty (no day-resolution tips)"

    if (cloud["rtt_mutations"] < 0).any():
        differences.append(
            f"root-to-tip: {int((cloud['rtt_mutations'] < 0).sum())} cell(s) have "
            f"negative divergence, which means the path accumulation ran the wrong "
            f"way; the scatter would still look like a cloud"
        )
    if (cloud["n_tips"] < 1).any():
        differences.append("root-to-tip: a cell holds fewer than 1 tip")

    n_cloud_tips = int(cloud["n_tips"].sum())
    residuals_path = os.path.join(output_dir, "mat", "node_dates_residuals.csv")
    if os.path.exists(residuals_path):
        n_residual_tips = sum(1 for _ in open(residuals_path)) - 1
        if n_cloud_tips > n_residual_tips:
            differences.append(
                f"root-to-tip: the cloud holds {n_cloud_tips} tips but only "
                f"{n_residual_tips} day-resolution tips exist; it can be smaller "
                f"(strains absent from this tree) but never larger"
            )

    # The report's own number, against a recomputation from the file it names.
    slope_reported = None
    if os.path.exists(report_path):
        report = pd.read_csv(report_path, sep="\t")
        row = report[report["check"] == "root_to_tip_regression"]
        if row.empty:
            differences.append(
                f"root-to-tip: no root_to_tip_regression row in {report_path}"
            )
        else:
            match = re.search(
                r"recomputed slope ([-\d.]+) subs/yr", str(row.iloc[0]["detail"])
            )
            if not match:
                differences.append(
                    f"root-to-tip: could not read the slope out of the report "
                    f"detail: {row.iloc[0]['detail']!r}"
                )
            else:
                slope_reported = float(match.group(1))

    if slope_reported is not None:
        qc = _load_node_date_qc_module()
        days = cloud["observed_date"].to_numpy(dtype="datetime64[D]").astype("float64")
        recomputed = qc.weighted_linregress(
            days,
            cloud["rtt_mutations"].to_numpy(dtype="float64"),
            cloud["n_tips"].to_numpy(dtype="float64"),
        )
        if recomputed["slope"] is None:
            differences.append(
                "root-to-tip: the cloud gives no slope but the report printed one"
            )
        else:
            per_year = recomputed["slope"] * qc.DAYS_PER_YEAR
            # The report rounds to 2 dp, so compare at that resolution.
            if abs(per_year - slope_reported) > 0.02:
                differences.append(
                    f"root-to-tip: the report says {slope_reported} subs/yr but "
                    f"{os.path.basename(cloud_path)} gives {per_year:.2f}. The "
                    f"figure regresses this same file, so the report and the "
                    f"figure are describing different regressions."
                )

    # 4. The figure's anchor for the fitted-clock line, against the root date the
    #    report states. Same script writes both, but they are read by different
    #    audiences -- the sidecar by the plot, the report by a person -- so a
    #    divergence would anchor the model-implied line somewhere nobody can see,
    #    and the line would still look like a line.
    anchor_path = os.path.join(output_dir, "mat", "node_dates_rtt_anchor.csv")
    anchor_note = "no anchor sidecar"
    if not os.path.exists(anchor_path):
        differences.append(
            f"root-to-tip: {anchor_path} is missing, so the fitted-clock line has "
            f"no root to start from"
        )
    else:
        anchor = pd.read_csv(anchor_path)
        if len(anchor) != 1:
            differences.append(
                f"root-to-tip: {anchor_path} holds {len(anchor)} rows, expected "
                f"exactly 1 (one root)"
            )
        else:
            anchor_date = pd.to_datetime(anchor["root_fitted_date"].iloc[0])
            anchor_note = f"anchor {anchor_date.date() if pd.notna(anchor_date) else 'blank'}"
            if os.path.exists(report_path):
                root_row = report[report["check"] == "root_date_plausible"]
                match = re.search(
                    r"dated (\d{4}-\d{2}-\d{2})", str(root_row.iloc[0]["detail"])
                ) if not root_row.empty else None
                if match and pd.notna(anchor_date):
                    if str(anchor_date.date()) != match.group(1):
                        differences.append(
                            f"root-to-tip: the anchor sidecar says the root is dated "
                            f"{anchor_date.date()} but check 7 reported "
                            f"{match.group(1)}; the figure's fitted-clock line is "
                            f"anchored somewhere the report does not describe"
                        )
            if pd.notna(anchor_date):
                # The root cannot postdate the earliest sample beneath it, so its
                # anchor sits at or left of the cloud. Inside the cloud would make
                # the model-implied line start partway up the y-axis at zero
                # divergence, which is not a rendering problem but a dating one.
                earliest = cloud["observed_date"].min()
                if anchor_date.normalize() > earliest:
                    differences.append(
                        f"root-to-tip: the root is dated {anchor_date.date()}, after "
                        f"the earliest tip in the cloud ({earliest.date()})"
                    )

    return differences, (
        f"root-to-tip cloud: {n_cloud_tips} tips in {len(cloud)} exact cells, "
        f"max divergence {int(cloud['rtt_mutations'].max())}, report slope "
        f"{slope_reported} subs/yr reproduced from the cloud, {anchor_note}"
    )


def _load_node_date_qc_module():
    """Import ``workflow/modules/node_date_qc.py`` by path (no package install)."""
    import importlib.util
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    module_path = os.path.join(
        os.path.dirname(scripts_dir), "modules", "node_date_qc.py"
    )
    if not os.path.exists(module_path):
        raise ImportError(f"Cannot find module 'node_date_qc' at {module_path}")
    spec = importlib.util.spec_from_file_location("node_date_qc", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["node_date_qc"] = module
    spec.loader.exec_module(module)
    return module


def _load_masking_module():
    """Import ``workflow/modules/masking.py`` by path (no package install).

    Mirrors the established ``load_workflow_module`` helper used by the pipeline
    scripts (e.g. ``aamut_fitness.py``): guards on a missing module file and
    registers the module in ``sys.modules``.
    """
    import importlib.util
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    module_path = os.path.join(
        os.path.dirname(scripts_dir), "modules", "masking.py"
    )
    if not os.path.exists(module_path):
        raise ImportError(f"Cannot find module 'masking' at {module_path}")
    spec = importlib.util.spec_from_file_location("masking", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["masking"] = module
    spec.loader.exec_module(module)
    return module


def _load_founder_base_fn(output_dir):
    """Build ``founder_nt(clade, site)`` from the run's reference + founders.

    Founder base = the reference nt at the site, overridden by the clade's own
    nucleotide substitutions (Neher primary, Roemer fallback) -- mirroring
    ``clade_founder_fasta.py``. Reads the FASTA by hand so no biopython dep is
    needed. Returns ``None`` if the required run outputs aren't present, and the
    lookup returns ``None`` for clades absent from both founder sources.
    """
    neher_path = os.path.join(output_dir, "clade_founders", "clade_founders_neher.json")
    roemer_path = os.path.join(output_dir, "clade_founders", "clade_founders_roemer.json")
    ref = _read_reference(output_dir)
    if ref is None or not os.path.exists(neher_path):
        return None

    with open(neher_path) as f:
        neher = json.load(f)
    roemer = {}
    if os.path.exists(roemer_path):
        with open(roemer_path) as f:
            roemer = json.load(f)

    def clade_subs(clade):
        info = neher.get(clade)
        if info is not None:
            return info.get("nuc") or info.get("mutations") or info.get("nuc_mutations") or []
        info = roemer.get(clade)
        if info is not None:
            return info.get("nucSubstitutions") or info.get("substitutions") or []
        return None

    def founder_nt(clade, site):
        subs = clade_subs(clade)
        if subs is None or site < 1 or site > len(ref):
            return None
        base = ref[site - 1]
        for mut in subs:
            try:
                pos = int(mut[1:-1])
            except (ValueError, IndexError):
                continue
            if pos == site:
                base = mut[-1]
        return base

    return founder_nt


def validate_usher_mutation_masking(output_dir, yaml_path="data/usher_masked_sites.yaml"):
    """Validate per-clade mutation masking in ``masking.compute_exclusions``.

    The pipeline's test clades (19A/19B/20A) carry no per-clade ``mutations:``
    masks, so ``counts_all.csv`` never exercises that path. This drives
    ``compute_exclusions`` directly against the committed UShER yaml with a
    witness reversion that must be:

      * masked in the clade that lists it (both ``exclude`` and
        ``masked_in_usher``),
      * NOT masked in a clade that lists neither the mutation nor its site
        (proves masking is clade-specific), and
      * masked at the mutation level, not the whole site: a *different*
        destination base at the same site, in the same clade, must stay
        unmasked (proves the mask is directional).

    The witness is chosen so the reversion is a genuine founder-relative grid
    row in *both* clades -- i.e. each clade's founder base at the site equals
    the mutation's "from" base -- so it is a real, countable mutation excluded
    in one clade but not the other, not merely a synthetic row. This requires
    the run's reference + founder outputs under ``output_dir``; if they are
    missing, or no founder-valid witness exists, the check fails (there is no
    weaker fallback).

    Returns ``(differences, summary)`` where ``summary`` describes the witness
    when the check passes (``None`` on failure).
    """
    if not os.path.exists(yaml_path):
        return [f"UShER masked-sites yaml not found at {yaml_path}"], None
    try:
        import pandas as pd
    except ImportError:
        return ["pandas unavailable; cannot validate mutation masking"], None

    masking = _load_masking_module()
    sites, muts = masking.load_usher_masks(yaml_path)

    if not any(muts.values()):
        return (
            [f"No per-clade `mutations:` masks found in {yaml_path}; "
             "mutation-level masking is unverified"],
            None,
        )

    def site_of(mut):
        return int(mut[1:-1])

    founder_nt = _load_founder_base_fn(output_dir)
    if founder_nt is None:
        return (
            [f"Reference/clade-founder outputs not found under {output_dir}; "
             "cannot verify mutation masking (need ref/ref.fa and "
             "clade_founders/clade_founders_neher.json)"],
            None,
        )
    all_clades = sorted(set(sites) | set(muts))

    def find_witness():
        """First founder-valid (clade_A, mut, site, clade_B), sorted for stability.

        clade_A lists the mutation but not its site; clade_B masks neither the
        mutation nor its site; and both clades' founder base at the site equals
        the mutation's "from" base, so the reversion is a real founder-relative
        grid row in each.
        """
        for ca in sorted(muts):
            for mut in sorted(muts[ca]):
                s = site_of(mut)
                frm = mut[0]
                if s in sites.get(ca, set()) or founder_nt(ca, s) != frm:
                    continue
                for cb in all_clades:
                    if cb == ca or mut in muts.get(cb, set()) or s in sites.get(cb, set()):
                        continue
                    if founder_nt(cb, s) != frm:
                        continue
                    return (ca, mut, s, cb)
        return None

    witness = find_witness()
    if witness is None:
        return (
            ["Could not find a founder-valid mutation-mask witness (a reversion "
             "that is a founder-relative grid row in both a clade that masks it "
             "and one that does not); mutation masking unverified"],
            None,
        )

    ca, mut, site, cb = witness
    from_nt, to_nt = mut[0], mut[-1]

    # A different destination base at the same site, avoiding any base that is
    # itself masked for clade_A at that site (in case of multiple reversions).
    masked_tos = {m[-1] for m in muts.get(ca, set()) if site_of(m) == site}
    other_to = next(
        (n for n in "ACGT" if n != from_nt and n not in masked_tos), None
    )

    rows = {
        "clade": [ca, cb],
        "nt_site": [site, site],
        "nt_mutation": [mut, mut],
    }
    if other_to is not None:
        other_mut = f"{from_nt}{site}{other_to}"
        rows["clade"].append(ca)
        rows["nt_site"].append(site)
        rows["nt_mutation"].append(other_mut)

    df = pd.DataFrame(rows)
    # reference_nt="N" disables the (separately-tested) founder->reference
    # reversion path so this check isolates the UShER mutation-mask mechanism.
    df["reference_nt"] = "N"
    exclude, masked = masking.compute_exclusions(
        df,
        sites_to_exclude=[],
        vcf_masked_sites=[],
        usher_clade_sites=sites,
        usher_clade_mutations=muts,
        site_include_range={"start": 266, "end": 29674},
    )
    exclude = [bool(x) for x in exclude]
    masked = [bool(x) for x in masked]

    differences = []
    # Row 0: masked in the clade that lists the mutation.
    if not (exclude[0] and masked[0]):
        differences.append(
            f"mutation mask: {mut} should be masked in {ca} "
            f"(exclude={exclude[0]}, masked_in_usher={masked[0]})"
        )
    # Row 1: not masked in a clade that lists neither the mutation nor its site.
    if exclude[1] or masked[1]:
        differences.append(
            f"mutation mask: {mut} should NOT be masked in {cb} "
            f"(exclude={exclude[1]}, masked_in_usher={masked[1]})"
        )
    # Row 2: a different destination base at the same site must stay unmasked
    # in clade_A (mutation masks are directional, not whole-site).
    if other_to is not None and (exclude[2] or masked[2]):
        differences.append(
            f"mutation mask: {other_mut} should NOT be masked in {ca}; "
            f"per-clade mutation masks must be directional, not whole-site "
            f"(exclude={exclude[2]}, masked_in_usher={masked[2]})"
        )

    summary = None if differences else (
        f"{mut} masked in {ca}, not in {cb} (founder-valid grid row in both)"
    )
    return differences, summary


def validate_usher_site_masking(yaml_path="data/usher_masked_sites.yaml"):
    """Validate per-clade UShER *site* masking in ``compute_exclusions``.

    Complements the mutation-level check: a whole-site mask must apply in the
    clade(s) that list it and not in a clade that does not. Founder identity is
    irrelevant here -- a site mask hits every mutation at the site -- so this
    drives ``compute_exclusions`` with an arbitrary mutation at the witness site
    and asserts the site is excluded (``exclude`` + ``masked_in_usher``) in one
    clade but not another.

    Returns ``(differences, summary)`` (``summary`` is ``None`` on failure).
    """
    if not os.path.exists(yaml_path):
        return [f"UShER masked-sites yaml not found at {yaml_path}"], None
    try:
        import pandas as pd
    except ImportError:
        return ["pandas unavailable; cannot validate site masking"], None

    masking = _load_masking_module()
    sites, muts = masking.load_usher_masks(yaml_path)

    if not any(sites.values()):
        return (
            [f"No per-clade `sites:` masks found in {yaml_path}; "
             "site-level masking is unverified"],
            None,
        )

    RANGE = {"start": 266, "end": 29674}
    all_clades = sorted(set(sites) | set(muts))

    # Witness: a site masked in clade_A but not in clade_B, kept inside the
    # include range so the range filter can't confound clade_B's result.
    witness = None
    for ca in sorted(sites):
        for s in sorted(sites[ca]):
            if not (RANGE["start"] <= s <= RANGE["end"]):
                continue
            for cb in all_clades:
                if cb == ca or s in sites.get(cb, set()):
                    continue
                witness = (ca, s, cb)
                break
            if witness:
                break
        if witness:
            break

    if witness is None:
        return (
            ["Could not find a site masked in one clade but not another; "
             "site masking unverified"],
            None,
        )

    ca, site, cb = witness
    # Any mutation at the site works; pick one that is not itself a mutation
    # mask in either clade so only the site mask can drive exclusion.
    mut = next(
        f"{a}{site}{b}"
        for a in "ACGT" for b in "ACGT"
        if a != b
        and f"{a}{site}{b}" not in muts.get(ca, set())
        and f"{a}{site}{b}" not in muts.get(cb, set())
    )

    df = pd.DataFrame(
        {"clade": [ca, cb], "nt_site": [site, site], "nt_mutation": [mut, mut]}
    )
    # reference_nt="N" disables the founder->reference reversion path so this
    # check isolates the whole-site UShER mask mechanism.
    df["reference_nt"] = "N"
    exclude, masked = masking.compute_exclusions(
        df,
        sites_to_exclude=[],
        vcf_masked_sites=[],
        usher_clade_sites=sites,
        usher_clade_mutations=muts,
        site_include_range=RANGE,
    )
    exclude = [bool(x) for x in exclude]
    masked = [bool(x) for x in masked]

    differences = []
    if not (exclude[0] and masked[0]):
        differences.append(
            f"site mask: site {site} should be masked in {ca} "
            f"(exclude={exclude[0]}, masked_in_usher={masked[0]})"
        )
    if exclude[1] or masked[1]:
        differences.append(
            f"site mask: site {site} should NOT be masked in {cb} "
            f"(exclude={exclude[1]}, masked_in_usher={masked[1]})"
        )

    summary = None if differences else f"site {site} masked in {ca}, not in {cb}"
    return differences, summary


def _read_reference(output_dir):
    """Read the run's reference FASTA (single record) as a plain string."""
    ref_path = os.path.join(output_dir, "ref", "ref.fa")
    if not os.path.exists(ref_path):
        return None
    seq = []
    with open(ref_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def validate_reversion_masking(output_dir):
    """Every founder->reference reversion in counts_all.csv must be excluded.

    A reversion is a grid row whose destination base equals the reference base
    at that site (``nt_mutation[-1] == reference_nt``). Checked against the real
    ``counts_all.csv`` (the test clades DO have clade-defining sites, e.g. 20A,
    so reversions are present). Also guards against vacuity (not every row is
    excluded, and reversions actually exist).
    """
    import pandas as pd

    counts_path = os.path.join(output_dir, "counts_all.csv")
    if not os.path.exists(counts_path):
        return [f"counts_all.csv not found under {output_dir}"], None
    ref = _read_reference(output_dir)
    if ref is None:
        return [f"ref/ref.fa not found under {output_dir}"], None

    df = pd.read_csv(counts_path)
    ref_nt = df["nt_site"].map(lambda s: ref[int(s) - 1])
    is_reversion = df["nt_mutation"].str[-1] == ref_nt
    n_rev = int(is_reversion.sum())
    if n_rev == 0:
        return ["no founder->reference reversion rows in counts_all.csv; "
                "reversion masking unverified"], None

    excluded = df["exclude"] == True  # noqa: E712
    leaked = df[is_reversion & ~excluded]
    differences = []
    if len(leaked):
        ex = leaked.iloc[0]
        differences.append(
            f"reversion masking: {len(leaked)} founder->reference reversion "
            f"row(s) are not excluded (e.g. clade {ex['clade']} {ex['nt_mutation']})"
        )
    if bool(excluded.all()):
        differences.append("reversion check is vacuous: every counts_all.csv row is excluded")

    summary = None if differences else (
        f"{n_rev} founder->reference reversion rows, all exclude=True"
    )
    return differences, summary


def _parse_codon_sites(series):
    """Flatten a codon_site column (handles ';'-joined and float values)."""
    out = set()
    for site in series.unique():
        if site is None:
            continue
        s = str(site)
        for part in s.split(";"):
            part = part.strip()
            if part.replace(".", "").replace("-", "").isdigit():
                out.add(int(float(part)))
    return out


def validate_spike_mask_all_clades(output_dir, masking_gene="S"):
    """The spike codon mask must equal the codon sites excluded in ALL clades.

    Re-derives the expected mask from counts_all.csv using the all-clades rule
    (a nt_site qualifies iff every row at it has exclude==True) restricted to the
    masking gene, and asserts the emitted ``{gene}_codon_sites_to_mask.txt`` matches.
    This pins the clade-agnostic behavior: a regression to per-clade masking would
    add clade-specific codon sites and fail here.
    """
    import pandas as pd

    counts_path = os.path.join(output_dir, "counts_all.csv")
    mask_path = os.path.join(output_dir, f"{masking_gene}_codon_sites_to_mask.txt")
    if not os.path.exists(counts_path):
        return [f"counts_all.csv not found under {output_dir}"], None
    if not os.path.exists(mask_path):
        return [f"{masking_gene}_codon_sites_to_mask.txt not found under {output_dir}"], None

    df = pd.read_csv(counts_path)
    site_all_excluded = (df["exclude"] == True).groupby(df["nt_site"]).transform("all")  # noqa: E712
    # Exact gene membership (mirrors get_problematic_sites.py, #114): split
    # compound cells (e.g. "S;ORF3a") on ";" so the mask can't bleed a gene
    # symbol into a longer one once this validator is used for other genes.
    gene_members = df["gene"].astype(str).str.split(";")
    in_gene = gene_members.apply(lambda members: masking_gene in [g.strip() for g in members])
    expected = _parse_codon_sites(df[site_all_excluded & in_gene]["codon_site"])

    with open(mask_path) as f:
        emitted = {int(x) for x in f.read().split() if x.strip().lstrip("-").isdigit()}

    differences = []
    if emitted != expected:
        differences.append(
            f"spike mask mismatch: emitted-only={sorted(emitted - expected)[:10]}, "
            f"expected-only={sorted(expected - emitted)[:10]}"
        )
    summary = None if differences else (
        f"{len(emitted)} codon sites, all excluded across every clade"
    )
    return differences, summary


_STOP_CODONS = {"TAA", "TAG", "TGA"}


def validate_orf1ab_founder(output_dir):
    """The ORF1ab gene-sliced founder must be a valid in-frame CDS.

    Confirms the two-block frameshift reconstruction (the GTF's ORF1ab CDS
    blocks share nt 13468, read twice): the 19A ORF1ab founder slice is 21288
    nt (7096 codons) and translates with no internal stop codon. Skips only when
    ORF1ab is not among ``pcp_genes`` (its founder subdir is absent), so
    Spike-only configs still pass; if the subdir exists but the file is missing
    that is a real failure, not a silent skip.

    Translation uses a hand-rolled codon table (stdlib only) since this script
    runs in an environment without biopython.
    """
    gene_dir = os.path.join(output_dir, "clade_founders", "ORF1ab")
    if not os.path.isdir(gene_dir):
        # ORF1ab not among pcp_genes for this run: nothing to validate.
        return [], None
    path = os.path.join(gene_dir, "19A.fa")
    if not os.path.exists(path):
        return [f"ORF1ab founder expected at {path} but missing "
                "(clade_founders/ORF1ab/ exists, so ORF1ab was configured)"], None

    with open(path) as f:
        seq = "".join(line.strip() for line in f if not line.startswith(">")).upper()

    differences = []
    if len(seq) % 3 != 0:
        differences.append(f"ORF1ab founder length {len(seq)} not divisible by 3")
    if len(seq) != 21288:
        differences.append(f"ORF1ab founder length {len(seq)} != expected 21288")
    stops = [
        (i // 3) + 1
        for i in range(0, len(seq) - len(seq) % 3, 3)
        if seq[i : i + 3] in _STOP_CODONS
    ]
    if stops:
        differences.append(
            f"ORF1ab founder translation has internal stop codon(s) at codon {stops[:5]}"
        )

    summary = None if differences else (
        f"ORF1ab founder in-frame: {len(seq)} nt = {len(seq) // 3} codons, no internal stop"
    )
    return differences, summary


def validate_orf1b_founder_slice(output_dir):
    """Every ORF1b founder must be exactly its clade's ORF1ab codons 4393-7096.

    ORF1b (nsp12-nsp16) is a PCP-only gene whose coordinates come from
    ``pcp_gene_cds_blocks`` rather than the reference GTF (#133). Its blocks,
    ``[[13442, 13468], [13468, 21552]]``, share nt 13468 exactly as ORF1ab's do,
    so a correct slice is byte-identical to the interior of the ORF1ab founder
    the pipeline builds independently. That makes this the cheapest available
    cross-check on the frameshift handling: the two founders are produced by the
    same code from different coordinates, and any frame error moves one and not
    the other.

    Checked for **every** clade, not just 19A, because a coordinate error would
    be identical across clades while a founder-construction error need not be.

    Skips when ORF1b is not among ``pcp_genes``. If ORF1b is present but ORF1ab
    is not, the identity cannot be checked; the frame checks still run and the
    summary says so rather than implying the comparison happened.
    """
    orf1b_dir = os.path.join(output_dir, "clade_founders", "ORF1b")
    if not os.path.isdir(orf1b_dir):
        # A visible SKIPPED, not a silent None: returning no summary makes a skip
        # indistinguishable from a pass in the "All validation checks passed"
        # block, which is how a check stops covering anything without anyone
        # noticing.
        #
        # ORF1b genuinely is optional in pcp_genes, so a skip can be correct --
        # but the *reason* is checked rather than asserted. If sibling gene
        # directories exist while ORF1b's does not, the run did request per-gene
        # founders and ORF1b's are missing, which is a failure and not a skip.
        founders_root = os.path.join(output_dir, "clade_founders")
        siblings = sorted(
            name for name in os.listdir(founders_root)
            if os.path.isdir(os.path.join(founders_root, name))
        ) if os.path.isdir(founders_root) else []
        if siblings:
            return [], (
                f"SKIPPED: no {orf1b_dir}; ORF1b is not among pcp_genes "
                f"(per-gene founders present for {', '.join(siblings)})"
            )
        return [], (
            f"SKIPPED: no per-gene founder directories under {founders_root} at "
            f"all, so ORF1b's absence is not evidence about pcp_genes"
        )
    orf1ab_dir = os.path.join(output_dir, "clade_founders", "ORF1ab")

    def read_seq(path):
        with open(path) as f:
            return "".join(
                line.strip() for line in f if not line.startswith(">")
            ).upper()

    differences = []
    fastas = sorted(f for f in os.listdir(orf1b_dir) if f.endswith(".fa"))
    if not fastas:
        return [f"clade_founders/ORF1b/ exists but holds no .fa files"], None

    n_compared = 0
    for name in fastas:
        clade = name[: -len(".fa")]
        seq = read_seq(os.path.join(orf1b_dir, name))
        if len(seq) != 8112:
            differences.append(
                f"ORF1b founder {clade} is {len(seq)} nt, expected 8112 "
                f"(2704 codons)"
            )
        stops = [
            (i // 3) + 1
            for i in range(0, len(seq) - len(seq) % 3, 3)
            if seq[i : i + 3] in _STOP_CODONS
        ]
        if stops:
            differences.append(
                f"ORF1b founder {clade} has internal stop codon(s) at "
                f"codon {stops[:5]}"
            )
        orf1ab_path = os.path.join(orf1ab_dir, name)
        if os.path.exists(orf1ab_path):
            orf1ab_seq = read_seq(orf1ab_path)
            expected = orf1ab_seq[13176:21288]
            if seq != expected:
                differences.append(
                    f"ORF1b founder {clade} != ORF1ab[13176:21288] "
                    f"({len(seq)} nt vs {len(expected)} nt); the frameshift "
                    f"reconstruction disagrees between the two coordinate sets"
                )
            else:
                n_compared += 1

    if differences:
        return differences, None
    if n_compared:
        summary = (
            f"ORF1b founders: {len(fastas)} clades, 8112 nt = 2704 codons, no "
            f"internal stop; {n_compared} byte-equal ORF1ab[13176:21288]"
        )
    else:
        summary = (
            f"ORF1b founders: {len(fastas)} clades in frame at 8112 nt, but "
            f"ORF1ab is not in pcp_genes so the slice identity was not checked"
        )
    return differences, summary


def validate_era_rate_tables(output_dir):
    """Check the per-era neutral rate tables against the run's own era summary.

    The era split (#137) fits one model per era and skips an era with no clades.
    Both halves of that are silent failures on their own: a missing
    ``master_table_{era}.csv`` looks the same whether the era was legitimately
    absent or was dropped by a wiring mistake, and a table that *is* present
    proves nothing about whether it was fitted on the right clades.

    ``glm/era_fit_summary.csv`` is the run's own record of which eras it fitted,
    so this reconciles it against what is actually on disk in both directions:
    every era marked ``fitted`` must have its two tables, and no era marked
    skipped may have them. A one-directional check would pass on a run that wrote
    every era's table and mislabelled the summary.

    On a single-era tree (every shipped test config, and the BA.2/JN.1/XFG trees)
    the era table must additionally be **identical** to the pooled one -- ISSUE-137
    criterion 5, single-era degeneracy, which is the guard that the split did not
    move the pre-existing fit. ``tests/test_era_rate_dispatch.py`` covers the
    two-era dispatch as a unit test, since no shipped config exercises it.
    """
    import csv

    differences = []
    summary_path = Path(output_dir) / 'glm' / 'era_fit_summary.csv'
    if not summary_path.exists():
        # The GLM stage may not have run at all (partial pipeline target).
        if (Path(output_dir) / 'glm' / 'master_table.csv').exists():
            differences.append(
                f"{summary_path} is missing although glm/master_table.csv exists, "
                "so which eras were fitted and which were skipped cannot be "
                "verified (#137)"
            )
        return differences, None

    with open(summary_path, newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        differences.append(f"{summary_path} has no rows")
        return differences, None

    fitted, skipped = [], []
    for row in rows:
        era = row['era']
        is_fitted = row['fitted'].strip().lower() == 'true'
        master = Path(output_dir) / 'glm' / f'master_table_{era}.csv'
        rates = Path(output_dir) / f'rates_table_{era}.csv'
        if is_fitted:
            fitted.append(era)
            for path in (master, rates):
                if not path.exists():
                    differences.append(
                        f"era {era} is marked fitted in {summary_path.name} but "
                        f"{path} does not exist"
                    )
            if int(row['n_clades']) == 0:
                differences.append(
                    f"era {era} is marked fitted with 0 clades, which cannot be "
                    "fitted"
                )
        else:
            skipped.append(era)
            for path in (master, rates):
                if path.exists():
                    differences.append(
                        f"era {era} is marked skipped in {summary_path.name} but "
                        f"{path} exists -- the summary and the output tree "
                        f"disagree about which eras were fitted"
                    )
            if int(row['n_clades']) != 0:
                differences.append(
                    f"era {era} is marked skipped but claims "
                    f"{row['n_clades']} clades"
                )

    # Single-era degeneracy (ISSUE-137 criterion 5): with one era, that era's
    # tables must equal the pooled ones exactly. This is the check that the
    # refactor did not move the v9 fit.
    if len(fitted) == 1 and not differences:
        era = fitted[0]
        for pooled_rel, era_rel in (
            (Path('glm') / 'master_table.csv', Path('glm') / f'master_table_{era}.csv'),
            (Path('rates_table.csv'), Path(f'rates_table_{era}.csv')),
        ):
            pooled = Path(output_dir) / pooled_rel
            per_era = Path(output_dir) / era_rel
            if not (pooled.exists() and per_era.exists()):
                continue
            if pooled.read_bytes() != per_era.read_bytes():
                differences.append(
                    f"this tree has only the {era} era, so {era_rel} must be "
                    f"identical to {pooled_rel} (ISSUE-137 criterion 5: "
                    f"single-era degeneracy). They differ, which means the era "
                    f"path is not the degenerate pooled path."
                )

    n_clades = sum(int(row['n_clades']) for row in rows)
    summary = (
        f"{len(fitted)} era(s) fitted ({', '.join(fitted) or 'none'}), "
        f"{len(skipped)} skipped ({', '.join(skipped) or 'none'}), "
        f"{n_clades} clades partitioned"
    )
    if len(fitted) == 1:
        summary += "; single-era tables byte-equal the pooled ones"
    return differences, summary


def validate_glm_synonymous_training(output_dir):
    """Check the GLM neutral model was trained on synonymous mutations only.

    ``glm/curated_counts.csv`` deliberately carries every mutation class -- the
    synonymous restriction is applied by ``modules/load.py::load_curated_muts``
    when ``master_table.py`` reads that file, not by ``curate_counts.py``. So
    the file's own composition proves nothing; what matters is how many rows the
    GLM actually consumed, which ``master_table.py`` logs.

    Non-synonymous and nonsense mutations are under purifying selection, so
    training on them biases the "neutral" rates downward and makes any
    actual-vs-expected fitness a comparison against average selection rather
    than against neutrality (#126).

    ``tests/test_glm_reference_equivalence.py`` covers the same contract as a
    unit test on a frozen fixture; this is the end-to-end counterpart that
    confirms the real pipeline run behaved.
    """
    import csv

    differences = []
    counts_path = Path(output_dir) / 'glm' / 'curated_counts.csv'
    log_path = Path(output_dir) / 'log' / 'master_table.log'

    # Skip only when the GLM stage never ran at all (partial pipeline target).
    # If the curated table exists, master_table must have run too, and a missing
    # log is a real problem: without it there is nothing to reconcile against,
    # so silently passing would make this check unenforceable.
    if not counts_path.exists():
        return differences, None
    if not log_path.exists():
        differences.append(
            f"{counts_path} exists but {log_path} does not, so the GLM's "
            "synonymous-only training cannot be verified (#126)"
        )
        return differences, None

    n_total = 0
    n_synonymous = 0
    with open(counts_path, newline='') as f:
        reader = csv.DictReader(f)
        if 'synonymous' not in (reader.fieldnames or []):
            differences.append(
                f"{counts_path} has no 'synonymous' column, so the GLM's "
                "synonymous-only filter cannot be verified (#126)"
            )
            return differences, None
        for row in reader:
            n_total += 1
            if row['synonymous'].strip().lower() == 'true':
                n_synonymous += 1

    with open(log_path) as f:
        log_text = f.read()
    match = re.search(r'Loaded (\d+) curated mutations for training', log_text)
    if not match:
        differences.append(
            f"{log_path} does not report how many mutations the GLM trained on; "
            "cannot verify the synonymous-only filter (#126)"
        )
        return differences, None
    n_trained = int(match.group(1))

    if n_trained != n_synonymous:
        differences.append(
            f"GLM trained on {n_trained} mutations but {counts_path.name} has "
            f"{n_synonymous} synonymous rows (of {n_total}). The neutral rate "
            "model must be trained on synonymous mutations only (#126)."
        )

    # If the curated table were somehow synonymous-only, the check above would
    # pass even with the filter removed. Flag that rather than pass vacuously.
    if n_synonymous == n_total:
        differences.append(
            f"{counts_path.name} is entirely synonymous ({n_total} rows), so it "
            "cannot distinguish a synonymous-only fit from a pooled one; the "
            "upstream curation or counting step has changed (#126)"
        )

    if differences:
        return differences, None
    return differences, (
        f"GLM trained on {n_trained}/{n_total} curated mutations, all synonymous"
    )


def validate_founder_codon_background(output_dir):
    """Check ``counts_all.csv`` holds exactly one row per (clade, nt_mutation).

    The grid enumerates every possible mutation from each clade founder, so this
    is a hard structural invariant, not a statistical one. It broke because the
    join carried no ``codon_change``: observations of the same nucleotide
    substitution on an already-mutated codon were relabelled onto the founder
    codon and kept as extra rows, corrupting both fitness estimators in opposite
    directions (#130). Nothing else could see it -- the phantom rows are
    byte-identical to the real ones in every emitted column.

    ``tests/test_founder_codon_background.py`` pins the filter itself on
    synthetic fixtures. This is the dataset-level counterpart, and it is the
    only one of the two that can speak about a *real* run.

    Note the test dataset has no off-founder-background observations at all, so
    a pass here confirms the invariant holds but does not exercise the filter;
    the reported summary says which case actually occurred.
    """
    import csv

    differences = []
    counts_path = Path(output_dir) / 'counts_all.csv'
    if not counts_path.exists():
        # A *difference*, not a skip: counts_all.csv is an unconditional output of
        # build_count_table for every tree, so there is no config under which it
        # is legitimately absent. Reporting SKIPPED here would pass the run on a
        # missing headline output.
        differences.append(
            f"{counts_path} does not exist, but build_count_table writes it for "
            f"every tree -- this is a missing output, not an inapplicable check"
        )
        return differences, None

    seen = set()
    duplicates = []
    n_rows = 0
    with open(counts_path, newline='') as f:
        reader = csv.DictReader(f)
        missing = [c for c in ('clade', 'nt_mutation') if c not in (reader.fieldnames or [])]
        if missing:
            differences.append(
                f"{counts_path.name} has no {missing} column(s), so the "
                "one-row-per-mutation invariant cannot be checked (#130)"
            )
            return differences, None
        for row in reader:
            n_rows += 1
            key = (row['clade'], row['nt_mutation'])
            if key in seen:
                if len(duplicates) < 5:
                    duplicates.append(key)
            else:
                seen.add(key)

    if duplicates:
        differences.append(
            f"{counts_path.name} has duplicated (clade, nt_mutation) rows -- "
            f"{n_rows} rows for {len(seen)} distinct mutations. Each mutation "
            "must appear exactly once per clade; extra rows mean observations "
            "from non-founder codon backgrounds were counted (#130). First "
            f"{len(duplicates)}: {duplicates}"
        )
        return differences, None

    # How many observations the filter actually removed, so a pass is not read
    # as "the filter was exercised" when the dataset had nothing to remove.
    log_path = Path(output_dir) / 'log' / 'build_count_table.log'
    dropped = 'log not found'
    if log_path.exists():
        match = re.search(
            r'dropped (\d+) rows mutated codon background', log_path.read_text()
        )
        dropped = f"{match.group(1)} dropped" if match else "none dropped"

    return differences, (
        f"Count table: {n_rows} rows, one per (clade, nt_mutation); "
        f"off-founder-background observations: {dropped}"
    )


def validate_cutoff_pcps(output_dir):
    """The date-cutoff PCP families reconcile with their roster, in both directions.

    Four properties, each of which fails silently without a check (#149).

    **The skip is real, not just recorded.** ``pcps/cutoff_summary.csv`` names
    every (cutoff, clade) pair with ``included`` true/false, and a false row must
    have *no* per-clade files while a true row must have all of them. Checking
    only the true rows would pass a run that also published a clade the
    checkpoint said to skip -- a subtree below ``min_clade_size``, or empty.
    Checking only the file tree would pass a run whose roster omitted a clade it
    built, which is the ambiguity ("no clade yet" vs "wiring mistake") the roster
    exists to resolve.

    **No cutoff PCP may involve a sample collected after its cutoff.** This is
    the property the whole feature is for, and nothing else here can see it: a
    date filter applied one comparison operator off, or against the wrong column,
    still yields plausible per-clade counts. Only leaf children can be checked --
    internal node names are not in the metadata -- but leaves are where the dates
    live.

    **No published sample is predicted well past its cutoff.** The
    reported-date check above cannot see the case the filter exists for: a metadata
    ``date`` that is simply wrong, which passes every date comparison because
    the comparison is against that same wrong value. So the published sample
    lists are re-checked against ``mat/node_dates.tsv`` using the bound the
    roster itself records, which is a recomputation rather than a re-reading of
    the count. Unlike the degenerate byte-equality below, this stays
    non-vacuous on production, where it is the only end-to-end statement that
    the filter did what the roster claims.

    **The degenerate arm is byte-identical to the full tree.** An arm where every
    clade's ``n_samples_upto`` equals its ``n_samples_full`` selected everything,
    so it must reproduce the full-tree aggregates exactly. That is detected from
    the summary rather than read from the config, so this function needs no
    knowledge of ``pcp_date_cutoffs``. Byte-equality (not row-equality) is the
    available check because ``pcpio.PCP_CSV_COMPRESSION`` pins ``mtime=0``
    (pattern 9); asserting only equal row counts would miss a reordering, and
    order is exactly what ``clades_w_adequate_counts_upto`` is written to
    preserve.

    If every arm turns out degenerate the result is reported as such rather than
    as a success, for #133's reason: "all cutoffs agreed" is true of a comparison
    that truncated nothing.
    """
    import pandas as pd

    summary_path = os.path.join(output_dir, "pcps", "cutoff_summary.csv")
    if not os.path.isfile(summary_path):
        return [], (
            f"SKIPPED: no {summary_path}; this config sets no pcp_date_cutoffs"
        )

    diffs = []
    summary = pd.read_csv(summary_path)
    expected_cols = {
        "cutoff", "clade", "output_subdir", "n_samples_full", "n_samples_upto",
        "n_predicted_late_dropped", "n_undated_kept", "max_predicted_days_past_cutoff",
        "adequate_sample_counts", "included",
    }
    missing_cols = expected_cols - set(summary.columns)
    if missing_cols:
        return [
            f"cutoff_summary.csv is missing column(s) {sorted(missing_cols)}"
        ], None
    if summary.empty:
        return [f"{summary_path} has no rows, but the file was requested"], None

    pcps_root = os.path.join(output_dir, "pcps")
    genes = sorted(
        name for name in os.listdir(pcps_root)
        if os.path.isdir(os.path.join(pcps_root, name))
    )
    if not genes:
        return [f"no per-gene PCP directories under {pcps_root}"], None

    # Directories on disk vs cutoffs in the roster, both ways. An orphan
    # ``upto_*`` directory is a stale output tree from a config that listed a
    # cutoff since removed, and its files would look as current as any other.
    roster_dirs = set(summary["output_subdir"])
    for gene in genes:
        on_disk = {
            name for name in os.listdir(os.path.join(pcps_root, gene))
            if name.startswith("upto_")
            and os.path.isdir(os.path.join(pcps_root, gene, name))
        }
        for orphan in sorted(on_disk - roster_dirs):
            diffs.append(
                f"pcps/{gene}/{orphan}/ exists but {orphan} is not in "
                f"cutoff_summary.csv (stale output tree from a removed cutoff?)"
            )
        for absent in sorted(roster_dirs - on_disk):
            diffs.append(
                f"cutoff_summary.csv names {absent} but pcps/{gene}/{absent}/ "
                f"does not exist"
            )

    # Sample dates, for the after-the-cutoff check below.
    meta_path = os.path.join(output_dir, "mat", "filtered_metadata.tsv.gz")
    strain_dates = None
    if os.path.isfile(meta_path):
        meta = pd.read_csv(
            meta_path, sep="\t", usecols=["strain", "date"], low_memory=False
        )
        strain_dates = dict(zip(meta["strain"], pd.to_datetime(meta["date"])))

    n_checked_pcps = 0
    n_late = 0
    for _, row in summary.iterrows():
        cutoff, clade, subdir = row["cutoff"], row["clade"], row["output_subdir"]
        included = bool(row["included"])
        if int(row["n_samples_upto"]) > int(row["n_samples_full"]):
            diffs.append(
                f"{cutoff}/{clade}: n_samples_upto "
                f"{row['n_samples_upto']} exceeds n_samples_full "
                f"{row['n_samples_full']}, which is arithmetically impossible"
            )

        samples_path = os.path.join(
            output_dir, "mat", "cutoffs", subdir, f"{clade}_samples.txt"
        )
        if included:
            if not os.path.isfile(samples_path):
                diffs.append(f"{cutoff}/{clade}: included but {samples_path} absent")
            else:
                with open(samples_path) as handle:
                    n_lines = sum(1 for line in handle if line.strip())
                if n_lines != int(row["n_samples_upto"]):
                    diffs.append(
                        f"{cutoff}/{clade}: sample list has {n_lines} strains but "
                        f"cutoff_summary.csv says {row['n_samples_upto']}"
                    )
        elif os.path.isfile(samples_path):
            # The same "recorded but not performed" test as the PCP check below,
            # one stage earlier. Checking only the PCP files leaves a whole class
            # of stale intermediate invisible: when the v12 production run first
            # mis-scored adequacy, the excluded clades' subtrees and sample lists
            # were built before the failure and survived the re-run that dropped
            # them, so the output tree described a fan-out that no longer existed.
            diffs.append(
                f"{cutoff}/{clade}: recorded included=False, yet "
                f"{samples_path} exists -- a stale intermediate from an "
                f"earlier fan-out, not a skip"
            )
        for gene in genes:
            pcp_path = os.path.join(pcps_root, gene, subdir, f"{clade}_pcps.csv.gz")
            if included and not os.path.isfile(pcp_path):
                diffs.append(f"{cutoff}/{clade}: included but {pcp_path} absent")
                continue
            if not included:
                if os.path.isfile(pcp_path):
                    diffs.append(
                        f"{cutoff}/{clade}: recorded included=False, yet "
                        f"{pcp_path} was published -- the skip was recorded but "
                        f"not performed"
                    )
                continue
            if strain_dates is None:
                continue
            pcps = pd.read_csv(pcp_path)
            bound = pd.Timestamp(cutoff)
            for name in pcps["child_name"]:
                date = strain_dates.get(name)
                if date is not None and date > bound:
                    n_late += 1
                    if n_late <= 5:
                        diffs.append(
                            f"{cutoff}/{clade}/{gene}: PCP child {name!r} was "
                            f"collected {date.date()}, after the cutoff"
                        )
            n_checked_pcps += len(pcps)

    # Independent recomputation of the predicted-date filter over what was actually
    # published. Reads the bound from the roster, so it cannot drift from the
    # config, and skips silently only when chronumental's output is absent --
    # which is reported, not swallowed.
    n_predicted_checked = 0
    n_predicted_late_published = 0
    node_dates_path = os.path.join(output_dir, "mat", "node_dates.tsv")
    predicted = None
    if os.path.isfile(node_dates_path):
        nd = pd.read_csv(
            node_dates_path,
            sep="\t",
            usecols=["strain", "predicted_date"],
            low_memory=False,
        )
        predicted = dict(
            zip(nd["strain"], pd.to_datetime(nd["predicted_date"]).dt.normalize())
        )
    for _, row in summary.iterrows():
        if not bool(row["included"]) or pd.isna(row["max_predicted_days_past_cutoff"]):
            continue
        if predicted is None:
            diffs.append(
                f"cutoff_summary.csv records a bound of "
                f"{row['max_predicted_days_past_cutoff']} days but "
                f"{node_dates_path} is absent, so nothing can check that the "
                f"published sets honour it"
            )
            break
        bound = pd.Timestamp(row["cutoff"]) + pd.Timedelta(
            days=int(row["max_predicted_days_past_cutoff"])
        )
        samples_path = os.path.join(
            output_dir, "mat", "cutoffs", row["output_subdir"],
            f"{row['clade']}_samples.txt",
        )
        if not os.path.isfile(samples_path):
            continue
        with open(samples_path) as handle:
            for line in handle:
                strain = line.strip()
                if not strain:
                    continue
                date = predicted.get(strain)
                if date is None:
                    continue
                n_predicted_checked += 1
                if date > bound:
                    n_predicted_late_published += 1
                    if n_predicted_late_published <= 5:
                        diffs.append(
                            f"{row['cutoff']}/{row['clade']}: published sample "
                            f"{strain!r} is predicted {date.date()}, more than "
                            f"{int(row['max_predicted_days_past_cutoff'])} days "
                            f"past the cutoff -- the filter did not remove it"
                        )

    n_late_total = int(summary["n_predicted_late_dropped"].sum())
    n_undated_total = int(summary["n_undated_kept"].sum())
    bounds = sorted(set(summary["max_predicted_days_past_cutoff"].dropna().astype(int)))
    if len(bounds) > 1:
        diffs.append(
            f"cutoff_summary.csv records more than one predicted-date bound "
            f"{bounds}; the config carries a single scalar, so the roster "
            f"mixes rows from different runs"
        )
    if n_late_total and not bounds:
        diffs.append(
            f"cutoff_summary.csv reports {n_late_total} predicted-late "
            f"drop(s) but no max_predicted_days_past_cutoff, so nothing records the bound "
            f"that produced the published sets"
        )
    grace_bound = bounds[0] if bounds else "none"

    # The degenerate arm(s), detected from the roster rather than the config.
    #
    # "Degenerate" means the *cutoff date* excluded nothing, which is not the
    # same as ``n_samples_upto == n_samples_full`` once the predicted-date
    # filter (#149) can also remove samples. Comparing the totals directly would
    # reclassify a degenerate arm as truncating the moment one predicted-late sample
    # was dropped, and the byte-equality check below -- the one check that says
    # the cutoff machinery left the ordinary path alone -- would simply stop
    # running. A check that silently stops running is pattern 10's vacuity.
    degenerate, truncating = [], []
    for cutoff, arm in summary.groupby("cutoff"):
        reached_cutoff = arm["n_samples_upto"] + arm["n_predicted_late_dropped"]
        if (reached_cutoff == arm["n_samples_full"]).all():
            degenerate.append(str(cutoff))
        else:
            truncating.append(str(cutoff))

    n_byte_equal = 0
    n_filter_confirmed = 0
    for cutoff in degenerate:
        subdir = f"upto_{cutoff}"
        arm = summary[summary["cutoff"].astype(str) == cutoff]
        n_dropped = int(arm["n_predicted_late_dropped"].sum())
        for gene in genes:
            for name in ("all_pcps.csv.gz", "all_pcps_dasm_ready.csv.gz"):
                full = os.path.join(pcps_root, gene, name)
                cut = os.path.join(pcps_root, gene, subdir, name)
                if not (os.path.isfile(full) and os.path.isfile(cut)):
                    continue
                with open(full, "rb") as a, open(cut, "rb") as b:
                    identical = a.read() == b.read()
                if n_dropped == 0:
                    if not identical:
                        diffs.append(
                            f"{cutoff} selected every sample and dropped none as "
                            f"predicted-late, so pcps/{gene}/{subdir}/{name} "
                            f"must byte-equal pcps/{gene}/{name}; it does not"
                        )
                    else:
                        n_byte_equal += 1
                elif identical:
                    # The converse direction, and the only end-to-end evidence
                    # that the filter reaches the extraction rather than merely
                    # being counted in the roster: a bound that removed samples
                    # cannot leave the outputs bit-for-bit unchanged.
                    diffs.append(
                        f"{cutoff} dropped {n_dropped} sample(s) as "
                        f"predicted-late, yet pcps/{gene}/{subdir}/{name} "
                        f"byte-equals pcps/{gene}/{name} -- the filter was "
                        f"recorded but did not reach the extraction"
                    )
                else:
                    n_filter_confirmed += 1

    # Truncating arms must actually shrink something, or the after-the-cutoff
    # check above ran over sets identical to the full tree's.
    for cutoff in truncating:
        arm = summary[summary["cutoff"].astype(str) == cutoff]
        lost_to_date = (
            arm["n_samples_upto"] + arm["n_predicted_late_dropped"]
            < arm["n_samples_full"]
        )
        if not lost_to_date.any():
            diffs.append(
                f"{cutoff} is classed as truncating but no clade lost a sample "
                f"to the cutoff date itself"
            )

    if not truncating:
        note = (
            f"VACUOUS: all {len(degenerate)} cutoff arm(s) selected every sample, "
            f"so nothing was truncated and the after-the-cutoff check compared "
            f"the full tree with itself"
        )
    else:
        skipped = summary[~summary["included"].astype(bool)]
        note = (
            f"{len(truncating)} truncating + {len(degenerate)} degenerate arm(s) "
            f"over {len(genes)} gene(s); {n_checked_pcps} PCP(s) checked, "
            f"{n_late} collected after their cutoff; {len(skipped)} "
            f"(cutoff, clade) pair(s) skipped and recorded; "
            f"{n_byte_equal} degenerate file(s) byte-equal the full tree; "
            f"{n_late_total} sample(s) dropped as predicted-late "
            f"(bound {grace_bound}), {n_undated_total} kept undated; "
            f"{n_predicted_checked} published sample(s) re-checked against "
            f"node_dates, {n_predicted_late_published} still late"
            + (
                f"; {n_filter_confirmed} degenerate file(s) differ as the "
                f"filter requires"
                if n_filter_confirmed
                else ""
            )
        )
    return diffs, note


def validate_subsample_arms(output_dir):
    """The subsampling arms reconcile with their roster, in both directions.

    Five properties, each of which fails silently without a check (#151).

    **The skip is real, not just recorded.** ``subsample_summary.csv`` names
    every (arm, clade) pair with ``included`` true/false, and a false row must
    have no arm subtree, no arm counts and no arm PCPs, while a true row must
    have all of them. Checking only the true rows would pass a run that also
    published a clade the sizing arithmetic said to thin away; checking only the
    file tree would pass a run whose roster omitted a clade it built -- the
    ambiguity ("thinned away" vs "wiring mistake") the roster exists to resolve.

    **An arm's samples are a subset of the clade's, of exactly the right size.**
    Set containment over the emitted lists, no duplicates, and the retained
    count against ``min(n_full, max_samples)``. This is check 1 of the issue,
    and it is the only statement that the sampler drew from the list it was
    pointed at.

    **An arm carries only clades strictly above its cap.** A clade at or
    below the cap is untouched, so an arm directory holding it would claim a
    thinning that did not happen and would turn one of the arm's own controls
    into an arm row. An arm that thins *nothing* (a cap at or above every
    clade's size) is legitimate and is recorded with an empty ``clade``: it must
    have no ``mat/clades/{arm}/`` directory at all, and being byte-identical to
    the full tree everywhere is the whole of its content.

    **A non-unit clade's per-clade files are the full tree's, byte-for-byte.**
    An arm re-aggregates every clade but re-counts only its unit, so the
    non-unit rows of ``mutation_counts/{arm}/aggregated.csv`` come from files
    that were never rewritten. Asserting the arm directory contains *only* unit
    clades is what makes that true by construction rather than by inspection --
    and it is what would catch a fan-out that routed an untouched clade through
    the sampler.

    **A degenerate arm is byte-identical to the full tree.** An arm whose every
    included clade retained all of its samples selected everything, so its
    aggregates must reproduce the full-tree ones exactly. Detected from the
    summary rather than read from the config, so this function needs no
    knowledge of ``clade_subsamples``. Byte-equality rather than row-equality
    because ``pcpio.PCP_CSV_COMPRESSION`` pins ``mtime=0`` (pattern 9), and
    because order is exactly what the draw's sort and ``arm_clades``' ordering
    are written to preserve -- equal row counts would miss a reordering.

    **The GLM was not refit.** No arm may have its own ``curated_counts*.csv``,
    ``master_table*.csv`` or ``rates_table*.csv``. The dependency chain
    ``annotate_counts -> curate_counts -> master_table -> rates_table`` makes
    refitting the default, so the cut has to be asserted rather than assumed:
    a refit arm would still produce well-formed tables, and its expected counts
    would silently be normalized against neutral rates fitted to the thinned
    data.

    If every arm turns out degenerate the result is reported as such rather than
    as a success, for #133's reason: "every arm agreed" is true of a comparison
    that thinned nothing.
    """
    import pandas as pd

    summary_path = os.path.join(output_dir, "subsample_summary.csv")
    if not os.path.isfile(summary_path):
        return [], (
            f"SKIPPED: no {summary_path}; this config sets no clade_subsamples"
        )

    diffs = []
    summary = pd.read_csv(summary_path)
    expected_cols = {
        "arm", "seed", "max_samples",
        "clade", "gene", "included",
        "n_samples_full", "n_samples_arm", "sample_fraction_realized",
        "n_nt_mutations_full", "n_nt_mutations_arm", "nt_mutation_retention",
    }
    missing_cols = expected_cols - set(summary.columns)
    if missing_cols:
        return [
            f"subsample_summary.csv is missing column(s) {sorted(missing_cols)}"
        ], None
    if summary.empty:
        return [f"{summary_path} has no rows, but the file was requested"], None

    # The counts-arm rows are the per-(arm, clade) roster; the gene rows repeat
    # each pair once per gene.
    roster = summary[summary["gene"].isna() | (summary["gene"] == "")]
    clades_root = os.path.join(output_dir, "mat", "clades")
    n_included = 0
    n_excluded = 0
    # Three kinds of arm, kept distinct because they support different claims:
    # a *degenerate* arm retained every sample of every clade it kept and so
    # must be byte-identical to the full tree; a *truncating* arm actually
    # thinned something and is the only kind that says the machinery does the
    # thing it exists for; an *excluding* arm dropped a clade below
    # ``min_clade_size`` and so differs from the full tree without having
    # thinned any surviving clade. Collapsing the last two into "not
    # degenerate" would let a run where nothing was ever thinned report success.
    degenerate_arms = []
    truncating_arms = set()
    excluding_arms = set()

    for arm, arm_rows in roster.groupby("arm"):
        # The unit is the set of clades the summary named for this arm. There is
        # no configured unit to cross-check it against: membership is derived
        # from the tree's clade sizes, which is why the summary records the
        # resolved list the fan-out actually used.
        unit = sorted(
            {str(c) for c in arm_rows["clade"] if str(c) and str(c) != "nan"}
        )
        arm_samples_dir = os.path.join(clades_root, arm)
        included_clades = []

        if not unit:
            # An arm that thins nothing: legitimate, and its whole content is
            # that it must not have drawn anything.
            if os.path.isdir(arm_samples_dir):
                diffs.append(
                    f"{arm}: summary records no thinned clade, but "
                    f"{arm_samples_dir} exists"
                )
            degenerate_arms.append(arm)
            continue

        for _, row in arm_rows.iterrows():
            clade = str(row["clade"])
            included = bool(row["included"])
            samples_path = os.path.join(arm_samples_dir, f"{clade}_samples.txt")
            exists = os.path.isfile(samples_path)

            if included and not exists:
                diffs.append(
                    f"{arm}/{clade}: summary says included but {samples_path} "
                    f"is missing"
                )
                continue
            if not included:
                n_excluded += 1
                excluding_arms.add(arm)
                if exists:
                    diffs.append(
                        f"{arm}/{clade}: summary says included=False but "
                        f"{samples_path} exists -- a stale intermediate would be "
                        f"read as a legitimate exclusion"
                    )
                # An excluded clade must also have no arm counts or PCPs.
                stale = os.path.join(
                    output_dir, "mutation_counts", arm, f"{clade}.csv"
                )
                if os.path.isfile(stale):
                    diffs.append(
                        f"{arm}/{clade}: summary says included=False but "
                        f"{stale} exists"
                    )
                continue

            n_included += 1
            included_clades.append(clade)
            full_path = os.path.join(clades_root, f"{clade}_samples.txt")
            with open(full_path) as f:
                full = [line.strip() for line in f if line.strip()]
            with open(samples_path) as f:
                drawn = [line.strip() for line in f if line.strip()]

            if not set(drawn) <= set(full):
                extra = sorted(set(drawn) - set(full))[:3]
                diffs.append(
                    f"{arm}/{clade}: drew {len(set(drawn) - set(full))} sample(s) "
                    f"not in the clade's list, e.g. {extra}"
                )
            if len(set(drawn)) != len(drawn):
                diffs.append(f"{arm}/{clade}: the drawn list has duplicates")

            cap = int(row["max_samples"])
            expected_n = min(len(full), cap)
            how = f"min({len(full)}, {cap})"
            # Strictly above the cap: a clade at or below it is untouched, so
            # its presence here means the arm claimed a clade it did not thin --
            # and that clade is no longer a control.
            if len(full) <= cap:
                diffs.append(
                    f"{arm}/{clade}: cap {cap} but the clade has only "
                    f"{len(full)} samples, so it is not thinned and must "
                    f"not appear in the arm"
                )
            if len(drawn) != expected_n:
                diffs.append(
                    f"{arm}/{clade}: drew {len(drawn)} samples, expected "
                    f"{how} = {expected_n}"
                )
            if int(row["n_samples_arm"]) != len(drawn):
                diffs.append(
                    f"{arm}/{clade}: summary records n_samples_arm="
                    f"{row['n_samples_arm']} but the file holds {len(drawn)}"
                )
            # Retention must run at or above the sample fraction: a branch
            # survives whenever any sample below it does. Strictly below would
            # mean mutations were lost that no collapse can explain.
            if row["n_nt_mutations_full"] and float(row["n_nt_mutations_full"]) > 0:
                retention = float(row["nt_mutation_retention"])
                realized = float(row["sample_fraction_realized"])
                if retention < realized - 1e-9:
                    diffs.append(
                        f"{arm}/{clade}: retained {retention:.4f} of the "
                        f"nucleotide mutations from {realized:.4f} of the "
                        f"samples; retention cannot fall below the sample "
                        f"fraction"
                    )
            if len(drawn) != len(full):
                truncating_arms.add(arm)

        # Only unit clades may have arm-local per-clade files. A non-unit clade
        # appearing here means the fan-out thinned a clade the config never
        # named, and its "unchanged" rows are not.
        if os.path.isdir(arm_samples_dir):
            present = sorted(
                name[: -len("_samples.txt")]
                for name in os.listdir(arm_samples_dir)
                if name.endswith("_samples.txt")
            )
            off_unit = [c for c in present if c not in unit]
            if off_unit:
                diffs.append(
                    f"{arm}: {arm_samples_dir} holds sample lists for non-unit "
                    f"clade(s) {off_unit}; only {unit} should be thinned"
                )

        # Degenerate only if it kept at least one clade, thinned none of them
        # and excluded none: an arm missing a clade from its grid cannot
        # reproduce the full tree's aggregates and must not be byte-compared.
        if (
            included_clades
            and arm not in truncating_arms
            and arm not in excluding_arms
        ):
            degenerate_arms.append(arm)

    # The GLM must not be refit per arm.
    for arm in sorted(roster["arm"].unique()):
        glm_arm_dir = os.path.join(output_dir, "glm", arm)
        if os.path.isdir(glm_arm_dir):
            refit = sorted(
                name for name in os.listdir(glm_arm_dir)
                if name.startswith(("curated_counts", "master_table", "rates_table"))
            )
            if refit:
                diffs.append(
                    f"{arm}: glm/{arm}/ holds {refit}, so the neutral model was "
                    f"refit on the thinned data. An arm must read the full "
                    f"tree's fit unchanged."
                )

    # Check 5: grid invariance. `build_count_table` enumerates the full mutation
    # grid from `clade_founder_nts.csv` and fills 0, so an arm's grid must differ
    # from the full tree's in the *count* columns and in nothing else. Comparing
    # every other column row-for-row catches a clade dropping out of the grid, a
    # founder changing, `site_include_range` applied differently -- and, because
    # `exclude` and `masked_in_usher` are among the compared columns, it is also
    # the check that shows the shared codon mask is a no-op rather than a leak.
    #
    # Restricted to the clades the arm kept: an *excluding* arm legitimately has
    # fewer clades, and comparing the whole file would report that intended
    # difference as a grid corruption. The clade set itself is already checked
    # above, so nothing is lost by keying on it here.
    COUNT_COLUMNS = [
        "actual_count",
        "count_terminal",
        "count_non_terminal",
        "mean_log_size",
    ]
    GRID_KEY = ["clade", "subset", "nt_site", "nt_mutation", "gene"]
    n_grids = 0
    full_grid_path = os.path.join(output_dir, "counts_all.csv")
    for arm in sorted(roster["arm"].unique()):
        arm_grid_path = os.path.join(output_dir, arm, "counts_all.csv")
        if not (os.path.isfile(full_grid_path) and os.path.isfile(arm_grid_path)):
            continue
        full_grid = pd.read_csv(full_grid_path)
        arm_grid = pd.read_csv(arm_grid_path)
        kept = set(arm_grid["clade"].unique())
        full_kept = full_grid[full_grid["clade"].isin(kept)]
        compare_cols = [c for c in full_grid.columns if c not in COUNT_COLUMNS]
        if list(arm_grid.columns) != list(full_grid.columns):
            diffs.append(
                f"{arm}: counts_all.csv columns differ from the full tree's "
                f"({list(arm_grid.columns)} vs {list(full_grid.columns)})"
            )
            continue
        if len(arm_grid) != len(full_kept):
            diffs.append(
                f"{arm}: counts_all.csv has {len(arm_grid)} rows over its "
                f"{len(kept)} clade(s), but the full tree has {len(full_kept)} "
                f"rows over the same clades. The grid is enumerated from the "
                f"founders and filled with 0, so it must not depend on "
                f"sampling depth."
            )
            continue
        n_grids += 1
        left = full_kept[compare_cols].sort_values(GRID_KEY).reset_index(drop=True)
        right = arm_grid[compare_cols].sort_values(GRID_KEY).reset_index(drop=True)
        unequal = [
            c for c in compare_cols
            if not left[c].equals(right[c])
            and not (left[c].isna() & right[c].isna() | (left[c] == right[c])).all()
        ]
        if unequal:
            first = unequal[0]
            bad = left.index[left[first] != right[first]]
            where = left.loc[bad[0], GRID_KEY].to_dict() if len(bad) else {}
            diffs.append(
                f"{arm}: counts_all.csv differs from the full tree's outside the "
                f"count columns, in {unequal}. First disagreement at {where}: "
                f"full={left.loc[bad[0], first] if len(bad) else '?'!r} vs "
                f"arm={right.loc[bad[0], first] if len(bad) else '?'!r}."
            )

    # Check 4: the per-clade normalization identity. `evol_time` is set so that
    # total expected synonymous count equals total observed synonymous count at
    # the training sites, and `genome_composition` counts *grid* rows -- which
    # check 5 just established are invariant to sampling depth. So
    # `K = sum(expected) / sum(actual)` over each clade's synonymous training
    # rows is a property of the grid and the shared rates alone, and must not
    # move with the cap. It is asserted per clade rather than
    # pooled: pooling lets one clade's K rise while another's falls.
    #
    # This is what makes an arm's fitness effects comparable to full depth at
    # all. A K that moved would mean the arm's expected counts were normalized
    # against a different total, so every `delta_fitness_ratio` would be shifted
    # by a factor no correlation against full depth could attribute.
    n_identities = 0
    full_exp_path = os.path.join(output_dir, "fitness", "expected_counts_by_clade.csv")

    def _clade_K(path):
        d = pd.read_csv(path)
        syn = d[(d["synonymous"] == True) & (d["exclude"] == False)]  # noqa: E712
        grouped = syn.groupby("clade")[["expected_count", "actual_count"]].sum()
        return {
            clade: row["expected_count"] / row["actual_count"]
            for clade, row in grouped.iterrows()
            if row["actual_count"] > 0
        }

    if os.path.isfile(full_exp_path):
        full_K = _clade_K(full_exp_path)
        for arm in sorted(roster["arm"].unique()):
            arm_exp_path = os.path.join(
                output_dir, "fitness", arm, "expected_counts_by_clade.csv"
            )
            if not os.path.isfile(arm_exp_path):
                continue
            arm_K = _clade_K(arm_exp_path)
            for clade, k_arm in sorted(arm_K.items()):
                if clade not in full_K:
                    diffs.append(
                        f"{arm}: clade {clade} has a normalization constant but "
                        f"the full tree does not, so the two cannot be compared"
                    )
                    continue
                k_full = full_K[clade]
                n_identities += 1
                if abs(k_arm - k_full) > 1e-10 * max(abs(k_full), 1.0):
                    diffs.append(
                        f"{arm}, clade {clade}: normalization constant "
                        f"sum(expected)/sum(actual) over synonymous training "
                        f"rows is {k_arm!r} against the full tree's {k_full!r} "
                        f"(relative difference "
                        f"{abs(k_arm - k_full) / abs(k_full):.3e}). The grid is "
                        f"invariant to sampling depth and the rates are shared, "
                        f"so this must not move."
                    )

    # A degenerate arm's aggregates must be byte-identical to the full tree's.
    n_compared = 0
    for arm in sorted(set(degenerate_arms)):
        pairs = [
            (os.path.join(output_dir, "mutation_counts", "aggregated.csv"),
             os.path.join(output_dir, "mutation_counts", arm, "aggregated.csv")),
            (os.path.join(output_dir, "counts_all.csv"),
             os.path.join(output_dir, arm, "counts_all.csv")),
            (os.path.join(output_dir, "aamut_fitness", "all_aamut_fitness.csv"),
             os.path.join(output_dir, "aamut_fitness", arm, "all_aamut_fitness.csv")),
        ]
        pcps_root = os.path.join(output_dir, "pcps")
        if os.path.isdir(pcps_root):
            for gene in sorted(os.listdir(pcps_root)):
                gene_dir = os.path.join(pcps_root, gene)
                if not os.path.isdir(gene_dir):
                    continue
                for name in ("all_pcps.csv.gz", "all_pcps_dasm_ready.csv.gz"):
                    pairs.append((
                        os.path.join(gene_dir, name),
                        os.path.join(gene_dir, arm, name),
                    ))
        for full_path, arm_path in pairs:
            if not (os.path.isfile(full_path) and os.path.isfile(arm_path)):
                continue
            n_compared += 1
            with open(full_path, "rb") as f:
                full_bytes = f.read()
            with open(arm_path, "rb") as f:
                arm_bytes = f.read()
            if full_bytes != arm_bytes:
                diffs.append(
                    f"degenerate arm {arm}: {os.path.relpath(arm_path, output_dir)} "
                    f"is not byte-identical to "
                    f"{os.path.relpath(full_path, output_dir)} "
                    f"({len(arm_bytes)} vs {len(full_bytes)} bytes)"
                )

    # correlation_summary.csv is the headline deliverable, and the two things a
    # run can check about it are that it exists for every arm and that its
    # coefficients rest on a population the reader can see. A coefficient
    # reported over an unstated n is the failure this table was designed against.
    corr_path = os.path.join(output_dir, "correlation_summary.csv")
    n_corr_thinned = 0
    corr_note = "absent"
    if not os.path.isfile(corr_path):
        diffs.append(
            f"correlation_summary.csv is missing from {output_dir}, so the arms "
            f"were built and never compared against full depth"
        )
    else:
        corr = pd.read_csv(corr_path)
        missing_arms = set(roster["arm"].unique()) - set(corr["arm"].unique())
        if missing_arms:
            diffs.append(
                f"correlation_summary.csv omits arm(s) {sorted(missing_arms)}, "
                f"which subsample_summary.csv lists"
            )
        # Every row must carry its n, including the empty ones: a row with a
        # coefficient and no n is unreadable, and a row silently dropped for
        # having an empty population is indistinguishable from one never asked.
        if corr["n_mutations"].isna().any():
            diffs.append(
                "correlation_summary.csv has row(s) with no n_mutations; a "
                "coefficient is not interpretable without the population it "
                "was computed on"
            )
        has_coef = corr["pearson_r"].notna()
        if (has_coef & (corr["n_mutations"].fillna(0) < 3)).any():
            bad = corr[has_coef & (corr["n_mutations"].fillna(0) < 3)]
            diffs.append(
                f"correlation_summary.csv reports {len(bad)} coefficient(s) "
                f"over fewer than 3 mutations, where Pearson r is undefined"
            )
        # A control clade reproducing full depth is the end-to-end evidence that
        # an arm did not reach a clade its config never named.
        control_moved = corr[(~corr["thinned"]) & (corr["rmse"].fillna(0) > 0)]
        if len(control_moved):
            diffs.append(
                f"correlation_summary.csv shows {len(control_moved)} control "
                f"row(s) (clades outside an arm's unit) whose fitness effects "
                f"differ from full depth: the arm reached a clade the config "
                f"never named"
            )
        thinned_rows = corr[corr["thinned"] & (corr["n_mutations"] > 0)]
        n_corr_thinned = len(thinned_rows)
        n_control = int((~corr["thinned"]).sum())
        maxcount = int(corr["max_actual_count_full"].max()) if len(corr) else 0
        if not n_corr_thinned:
            corr_note = (
                f"{len(corr)} row(s) but VACUOUS -- no thinned (arm, clade) "
                f"pair retained a scored mutation to compare"
            )
        else:
            corr_note = (
                f"{len(corr)} row(s), {n_corr_thinned} populated thinned "
                f"row(s), {n_control} control row(s) all matching full depth"
            )
            if maxcount <= 1:
                corr_note += (
                    f"; NOT A NOISE MEASUREMENT -- max full-depth actual_count "
                    f"is {maxcount}, so the coefficients are arithmetic"
                )

    n_arms = int(roster["arm"].nunique())
    if not truncating_arms:
        note = (
            f"VACUOUS: none of the {n_arms} arm(s) thinned a surviving clade, so "
            f"nothing here compares a subsampled arm against full depth "
            f"({len(degenerate_arms)} degenerate, {len(excluding_arms)} "
            f"excluding a clade)"
        )
    elif not degenerate_arms:
        note = (
            f"{n_arms} arm(s): {len(truncating_arms)} truncating, but none "
            f"degenerate -- no arm pins the byte-identity of the ordinary path "
            f"({n_included} included and {n_excluded} excluded pair(s))"
        )
    else:
        note = (
            f"{n_arms} arm(s): {len(truncating_arms)} truncating, "
            f"{len(degenerate_arms)} degenerate "
            f"({n_compared} file(s) byte-compared), {len(excluding_arms)} "
            f"excluding a clade; {n_included} included and {n_excluded} "
            f"excluded (arm, clade) pair(s); {n_grids} grid(s) invariant "
            f"outside the count columns; {n_identities} per-clade "
            f"normalization constant(s) unmoved; correlations: {corr_note}"
        )
    return diffs, note


def validate_censor_table(output_dir, yaml_path="data/usher_masked_sites.yaml"):
    """Every PCP gene's censor table is well formed and its roster is complete.

    The censor table (#144) carries the two *clade-specific* UShER reasons for
    DASM training. On the test data it is **empty for every gene** -- ``max_clades:
    3`` gives 19A/19B/20A and the yaml's earliest clade is 20B -- so the row-level
    properties are covered by ``tests/test_censor_table.py`` and cannot be
    asserted here. What a run *can* assert, and what is non-vacuous even when the
    table is empty, is the roster: a clade with nothing censored contributes zero
    rows, so the sidecar is the only thing that distinguishes it from a mistyped
    clade name. It must name every clade the run processed -- checked as an
    equality against ``clade_founders/{gene}/*.fa`` rather than as a subset of
    clade_counts.csv, since the latter is a superset of what any run processes
    and so cannot see a roster that drops a clade.

    Also checks that the two masking outputs are both present. They carry
    different exclusion reasons at different granularity, and a training run must
    apply **both** -- the censor table is not a replacement for
    ``{gene}_codon_sites_to_mask.txt``. Note the two are deliberately *not*
    asserted disjoint: a site can be both clade-masked in the yaml and in the
    ProblematicSites VCF (32 such pairs on the v11 gisaid-and-public tree), so
    disjointness is false on real data.
    """
    import pandas as pd

    pcps_root = os.path.join(output_dir, "pcps")
    if not os.path.isdir(pcps_root):
        return [], f"SKIPPED: no {pcps_root}, so no censor tables to validate"

    genes = sorted(
        name for name in os.listdir(pcps_root)
        if os.path.isdir(os.path.join(pcps_root, name))
    )
    if not genes:
        return [], f"SKIPPED: no per-gene PCP directories under {pcps_root}"

    # Two clade sets, checked in two different directions, because they fail
    # differently.
    #
    # `adequate_clades` is every clade clade_counts.csv calls adequate. It is a
    # *superset* of what a run processes -- `clades_w_adequate_counts` then
    # applies `clades_to_exclude` and truncates to `max_clades`, and both are set
    # in shipped configs (config.yaml excludes 21A; config_test_data.yaml caps at
    # 3) -- so it can only be used to catch a roster naming a clade that is not a
    # clade at all. Asserting the roster equals it would fail every correct run.
    #
    # `clade_founders/{gene}/*.fa` is the set the run actually processed for this
    # gene: the Snakefile materializes one founder slice per
    # `clades_w_adequate_counts` clade, which is the same list `censor_table.py`
    # writes as its roster. That admits an *equality* check, and equality is the
    # point -- the roster's whole job is to make "clade processed, nothing
    # censored" distinguishable from "clade never processed", which a
    # subset-of-a-superset test cannot see. Without it a roster silently missing
    # a processed clade passes, which is the vacuity this function exists to
    # prevent.
    clade_counts = os.path.join(output_dir, "mat", "clade_counts.csv")
    adequate_clades = None
    if os.path.isfile(clade_counts):
        counts = pd.read_csv(clade_counts)
        clade_col = next(
            (c for c in counts.columns if c.startswith("Nextstrain_clade")), None
        )
        if clade_col is not None:
            adequate_clades = set(
                counts[counts["adequate_sample_counts"] == True][clade_col].dropna()
            )

    def processed_clades_for(gene):
        """The clades this run built founder slices for, or None if unavailable."""
        gene_founders = os.path.join(output_dir, "clade_founders", gene)
        if not os.path.isdir(gene_founders):
            return None
        return {
            name[: -len(".fa")]
            for name in os.listdir(gene_founders)
            if name.endswith(".fa")
        }

    differences = []
    total_rows = 0
    roster_size = 0
    for gene in genes:
        table_path = os.path.join(pcps_root, gene, "censored_mutations.csv")
        roster_path = os.path.join(pcps_root, gene, "censored_mutations.clades.txt")
        codon_mask = os.path.join(output_dir, f"{gene}_codon_sites_to_mask.txt")

        if not os.path.isfile(table_path):
            differences.append(f"censor table missing for {gene}: {table_path}")
            continue
        if not os.path.isfile(roster_path):
            differences.append(f"censor clade roster missing for {gene}: {roster_path}")
            continue
        if not os.path.isfile(codon_mask):
            differences.append(
                f"{gene} has a censor table but no {codon_mask}. A training run "
                f"must apply both -- the censor table carries only the two "
                f"clade-specific UShER reasons, not the VCF/config/range ones."
            )

        # keep_default_na=False is the file's contract: a whole-site mask writes
        # empty strings, and pandas' defaults would turn them into NaN so that a
        # `from_nt == ""` test silently matches nothing.
        table = pd.read_csv(table_path, keep_default_na=False)
        expected_columns = ["clade", "nt_site", "nt_index", "from_nt", "to_nt"]
        if list(table.columns) != expected_columns:
            differences.append(
                f"{gene} censor table columns are {list(table.columns)}, expected "
                f"{expected_columns}"
            )
            continue
        total_rows += len(table)

        blank_from = table["from_nt"] == ""
        blank_to = table["to_nt"] == ""
        if not (blank_from == blank_to).all():
            differences.append(
                f"{gene} censor table has rows with exactly one of from_nt/to_nt "
                f"populated; a row must have both or neither"
            )

        roster = [line for line in open(roster_path).read().splitlines() if line]
        roster_size = max(roster_size, len(roster))
        if roster != sorted(roster):
            differences.append(f"{gene} censor clade roster is not sorted")
        if len(roster) != len(set(roster)):
            differences.append(f"{gene} censor clade roster has duplicates")
        if adequate_clades is not None and set(roster) - adequate_clades:
            differences.append(
                f"{gene} censor clade roster names clades absent from "
                f"clade_counts.csv: {sorted(set(roster) - adequate_clades)}"
            )

        processed = processed_clades_for(gene)
        if processed is None:
            differences.append(
                f"{gene} has a censor table but no clade_founders/{gene}/ "
                f"directory, so its roster cannot be checked for completeness. "
                f"The rule declares those founders as inputs, so they must exist "
                f"wherever the table does."
            )
        elif set(roster) != processed:
            missing = sorted(processed - set(roster))
            extra = sorted(set(roster) - processed)
            differences.append(
                f"{gene} censor clade roster does not match the clades this run "
                f"processed (clade_founders/{gene}/*.fa): missing {missing}, "
                f"unexpected {extra}. A clade omitted from the roster is "
                f"indistinguishable to a consumer from one that was never "
                f"processed, so 'nothing censored here' would read as 'unknown "
                f"clade'."
            )
        missing_from_roster = set(table["clade"]) - set(roster)
        if missing_from_roster:
            differences.append(
                f"{gene} censor table has rows for clades absent from its roster: "
                f"{sorted(missing_from_roster)}. The roster is what lets a "
                f"consumer tell an empty clade from a mistyped one, so it must "
                f"be a superset of the rows."
            )

    # The roster size is reported alongside the row count so that "0 rows" reads
    # as "checked N clades and none had a mask" rather than as "checked nothing".
    return differences, (
        f"Censor tables: {len(genes)} gene(s), {total_rows} row(s) total, "
        f"{roster_size} clade(s) per roster"
        + (" (empty on the test clades, as expected)" if total_rows == 0 else "")
    )


def default_output_dir():
    """Best-guess test output dir when none is given.

    The Snakefile appends the ``VERSION`` to ``output_dir`` (e.g.
    ``results/test_data`` -> ``results/test_data_v3``), so the bare
    ``results/test_data`` never exists. Prefer the versioned directory, then
    any ``results/test_data*``, falling back to the unversioned name.
    """
    import glob
    versioned = None
    try:
        with open('VERSION') as vf:
            versioned = f"results/test_data_{vf.read().strip()}"
        if os.path.isdir(versioned):
            return versioned
    except OSError:
        pass
    candidates = sorted(glob.glob('results/test_data*'))
    fallback = candidates[-1] if candidates else 'results/test_data'
    # Announce the fallback loudly (#88). Immediately after a version bump the
    # VERSION directory does not exist yet, so this silently validated the
    # *previous* version's outputs and reported success -- the failure #88 is
    # about. The fallback is kept (validating an older run on purpose is
    # legitimate) but it can no longer happen quietly. Pass an explicit path or
    # OUTPUT_DIR to bypass this entirely.
    if versioned and fallback != versioned:
        print(
            f"WARNING: {versioned} (from VERSION) does not exist; falling back to "
            f"{fallback}. These are DIFFERENT pipeline versions -- a pass here says "
            f"nothing about {versioned}. Re-run the pipeline, or pass an explicit "
            f"output directory.",
            file=sys.stderr,
        )
    return fallback


def main():
    """Run validation."""
    import os
    import subprocess

    # Determine output directory from args, env var, or default
    output_dir = (
        sys.argv[1] if len(sys.argv) > 1
        else os.environ.get('OUTPUT_DIR') or default_output_dir()
    )
    # Always name the directory being validated: "all checks passed" is
    # meaningless without knowing which version produced the outputs (#88).
    print(f"Validating {output_dir}")

    # Generate current summary
    result = subprocess.run(
        ['python', 'workflow/scripts/test_summary.py', output_dir],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"Failed to generate test summary: {result.stderr}")
        sys.exit(1)
    
    actual = json.loads(result.stdout)
    
    # Load expected summary
    expected_path = Path('data/test_data/expected_summary.json')
    if not expected_path.exists():
        print(f"Expected summary not found at {expected_path}")
        print("Current summary:")
        print(json.dumps(actual, indent=2))
        sys.exit(1)
    
    expected = load_json(expected_path)
    
    # Compare
    differences = compare_summaries(actual, expected)

    # Independently validate per-clade UShER masking against the committed yaml
    # (the test clades carry no UShER masks): both whole-site masks and the
    # optional directional mutation masks, each applied in one clade but not
    # another.
    site_mask_diffs, site_mask_summary = validate_usher_site_masking()
    mut_mask_diffs, mut_mask_summary = validate_usher_mutation_masking(output_dir)
    reversion_diffs, reversion_summary = validate_reversion_masking(output_dir)
    spike_mask_diffs, spike_mask_summary = validate_spike_mask_all_clades(output_dir)
    orf1ab_diffs, orf1ab_summary = validate_orf1ab_founder(output_dir)
    orf1b_diffs, orf1b_summary = validate_orf1b_founder_slice(output_dir)
    glm_syn_diffs, glm_syn_summary = validate_glm_synonymous_training(output_dir)
    codon_bg_diffs, codon_bg_summary = validate_founder_codon_background(output_dir)
    era_diffs, era_summary = validate_era_rate_tables(output_dir)
    censor_diffs, censor_summary = validate_censor_table(output_dir)
    rtt_diffs, rtt_summary = validate_root_to_tip_cloud(output_dir)
    cutoff_diffs, cutoff_summary = validate_cutoff_pcps(output_dir)
    arm_diffs, arm_summary = validate_subsample_arms(output_dir)
    differences += (
        site_mask_diffs + mut_mask_diffs + reversion_diffs + spike_mask_diffs
        + orf1ab_diffs + orf1b_diffs + glm_syn_diffs + codon_bg_diffs + era_diffs
        + censor_diffs + rtt_diffs + cutoff_diffs + arm_diffs
    )

    if differences:
        print("VALIDATION FAILED!")
        print("\nDifferences found:")
        for diff in differences:
            print(f"  - {diff}")
        print("\nActual summary:")
        print(json.dumps(actual, indent=2))
        sys.exit(1)
    else:
        print("✓ All validation checks passed!")
        print(f"  - Total PCPs: {actual['pcp_summary']['total_pcps']}")
        print(f"  - Total clades: {actual['clade_counts_summary']['total_clades']}")
        print(f"  - Total sequences: {actual['clade_counts_summary']['total_sequences']}")
        if actual.get('rates_table_summary'):
            print(f"  - Rates table rows: {actual['rates_table_summary']['num_rows']}")
        if actual.get('clade_founder_nts_summary'):
            print(f"  - Clade founder nucleotides: {actual['clade_founder_nts_summary']['num_rows']}")
        if actual.get('aamut_fitness_summary'):
            print(f"  - AA mutation fitness rows: {actual['aamut_fitness_summary']['num_rows']}")
            print(f"  - AA mutation fitness clade files: {actual['aamut_fitness_summary']['num_clade_files']}")
        if site_mask_summary:
            print(f"  - UShER per-clade site mask: {site_mask_summary}")
        if mut_mask_summary:
            print(f"  - UShER per-clade mutation mask: {mut_mask_summary}")
        if reversion_summary:
            print(f"  - Founder->reference reversions: {reversion_summary}")
        if spike_mask_summary:
            print(f"  - Spike codon mask (all-clades): {spike_mask_summary}")
        if orf1ab_summary:
            print(f"  - {orf1ab_summary}")
        if orf1b_summary:
            print(f"  - {orf1b_summary}")
        if glm_syn_summary:
            print(f"  - {glm_syn_summary}")
        if codon_bg_summary:
            print(f"  - {codon_bg_summary}")
        if era_summary:
            print(f"  - Era rate tables: {era_summary}")
        if censor_summary:
            print(f"  - {censor_summary}")
        if rtt_summary:
            print(f"  - {rtt_summary}")
        if cutoff_summary:
            print(f"  - Date-cutoff PCPs: {cutoff_summary}")
        if arm_summary:
            print(f"  - Clade-subsampling arms: {arm_summary}")
        sys.exit(0)


if __name__ == "__main__":
    main()