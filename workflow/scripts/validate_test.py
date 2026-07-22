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
            # to filter_viral_pcps.sh).
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
            for key in ['expected_count_sum', 'delta_fitness_pseudo_sum',
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
    in_gene = df["gene"].astype(str).str.contains(masking_gene, na=False)
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
    nt (7096 codons) and translates with no internal stop codon. Skips silently
    when ORF1ab is not among ``pcp_genes`` (so Spike-only configs still pass).

    Translation uses a hand-rolled codon table (stdlib only) since this script
    runs in an environment without biopython.
    """
    path = os.path.join(output_dir, "clade_founders", "ORF1ab", "19A.fa")
    if not os.path.exists(path):
        return [], None

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


def default_output_dir():
    """Best-guess test output dir when none is given.

    The Snakefile appends the ``VERSION`` to ``output_dir`` (e.g.
    ``results/test_data`` -> ``results/test_data_v3``), so the bare
    ``results/test_data`` never exists. Prefer the versioned directory, then
    any ``results/test_data*``, falling back to the unversioned name.
    """
    import glob
    try:
        with open('VERSION') as vf:
            versioned = f"results/test_data_{vf.read().strip()}"
        if os.path.isdir(versioned):
            return versioned
    except OSError:
        pass
    candidates = sorted(glob.glob('results/test_data*'))
    return candidates[-1] if candidates else 'results/test_data'


def main():
    """Run validation."""
    import os
    import subprocess

    # Determine output directory from args, env var, or default
    output_dir = (
        sys.argv[1] if len(sys.argv) > 1
        else os.environ.get('OUTPUT_DIR') or default_output_dir()
    )

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
    differences += (
        site_mask_diffs + mut_mask_diffs + reversion_diffs + spike_mask_diffs + orf1ab_diffs
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
        sys.exit(0)


if __name__ == "__main__":
    main()