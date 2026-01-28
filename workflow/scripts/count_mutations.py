#!/usr/bin/env python3
"""Count mutations from matUtils translated mutations."""

import Bio.SeqIO
import numpy as np
import pandas as pd
import yaml

# Get variables from snakemake
max_nt_mutations = snakemake.params.max_nt_mutations
max_reversions_to_ref = snakemake.params.max_reversions_to_ref
max_reversions_to_founder = snakemake.params.max_reversions_to_clade_founder
input_tsv = snakemake.input.tsv
input_nt_mut_csv = snakemake.input.nt_mut_csv
ref_fasta = snakemake.input.ref_fasta
usher_masked_sites_yaml = snakemake.input.usher_masked_sites
site_mask_csv = snakemake.input.site_mask
clade_founder_fasta = snakemake.input.clade_founder_fasta
clade = snakemake.wildcards.clade
sites_to_exclude_param = snakemake.params.sites_to_exclude
site_include_range = snakemake.params.site_include_range
exclude_ref_to_founder_muts = snakemake.params.exclude_ref_to_founder_muts
ref_to_founder_muts_csv = snakemake.input.ref_to_founder_muts
filtered_branches_csv = snakemake.input.filtered_branches
output_csv = snakemake.output.csv
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None


def log(msg):
    """Write to log file if available, otherwise print."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


log(f"Counting mutations for clade: {clade}")

# Get reference and founder sequence
ref = str(Bio.SeqIO.read(ref_fasta, "fasta").seq)
founder = str(Bio.SeqIO.read(clade_founder_fasta, "fasta").seq)
log(f"Reference length: {len(ref)}, Founder length: {len(founder)}")

# Build sites to exclude set
sites_to_exclude = set(sites_to_exclude_param) if sites_to_exclude_param else set()
masked_sites = set(pd.read_csv(site_mask_csv)["site"])
sites_to_exclude = sites_to_exclude.union(masked_sites)
log(f"Sites to exclude: {len(sites_to_exclude)}")

if exclude_ref_to_founder_muts:
    muts_to_exclude = set(pd.read_csv(ref_to_founder_muts_csv)["mutation"])
    log(f"Ref-to-founder mutations to exclude: {len(muts_to_exclude)}")
else:
    muts_to_exclude = set()

# Apply UShER masks for this clade
with open(usher_masked_sites_yaml) as f:
    usher_masked_sites = yaml.safe_load(f)

for mask, mask_dict in usher_masked_sites.items():
    if clade in mask_dict.get("clades", []):
        sites_to_exclude = sites_to_exclude.union(mask_dict.get("sites", []))
        log(
            f"Applied UShER mask '{mask}' with {len(mask_dict.get('sites', []))} sites"
        )

# Load filtered branches (branches with invalid PCPs to exclude)
filtered_branches_df = pd.read_csv(filtered_branches_csv)
filtered_branch_ids = set(filtered_branches_df["child_name"])
log(f"Filtered branches (invalid PCPs): {len(filtered_branch_ids)}")

# Expected final columns
final_cols = [
    "protein",
    "aa_mutation",
    "nt_mutation",
    "codon_change",
    "synonymous",
    "noncoding",
    "nt_site",
    "reference_nt",
    "clade_founder_nt",
    "exclude",
    "count",
    "count_terminal",
    "count_non_terminal",
    "mean_log_size",
]


def safe_mutation_count(mutations):
    """Count unique mutations, handling empty or malformed input."""
    if not mutations or mutations == [""]:
        return 0
    return len(set(mutations))


def count_reversions(mutations, reference_seq):
    """Count reversions to reference, with error handling."""
    if not mutations or mutations == [""]:
        return 0
    try:
        return sum(
            m[-1] == reference_seq[int(m[1:-1]) - 1]
            for m in set(mutations)
        )
    except (ValueError, IndexError) as e:
        log(f"Warning: Error counting reversions for mutations {mutations}: {e}")
        return 0


def get_noncoding(row):
    """Get just the noncoding mutations in a DataFrame row."""
    coding = row["nt_mutations"]
    both = row["nt_mutations_coding_and_noncoding"]
    assert set(coding).issubset(both)
    return sorted(set(both) - set(coding))


# Read input data
translated_mat = pd.read_csv(input_tsv, sep="\t")
nt_muts = pd.read_csv(input_nt_mut_csv).rename(
    columns={"nt_mutations": "nt_mutations_coding_and_noncoding"}
)

log(f"Translated MAT: {len(translated_mat)} rows")
log(f"NT mutations: {len(nt_muts)} rows")

if len(translated_mat) == 0 or len(nt_muts) == 0:
    log("No data, creating empty output")
    mutation_counts = pd.DataFrame(columns=final_cols)
else:
    # Filter nodes by mutation criteria
    filtered_nodes = (
        translated_mat.merge(nt_muts, on="node_id", validate="one_to_one")
        .query("node_id not in @filtered_branch_ids")  # Exclude branches with invalid PCPs
        .query("not nt_mutations.str.contains(',')")
        .assign(
            nt_mutations=lambda x: x["nt_mutations"].str.split(";"),
            nt_mutations_coding_and_noncoding=lambda x: x[
                "nt_mutations_coding_and_noncoding"
            ].str.split(";"),
            codon_changes=lambda x: x["codon_changes"].str.split(";"),
            aa_mutations=lambda x: x["aa_mutations"].str.split(";"),
            noncoding_nt_mutations=lambda x: x.apply(get_noncoding, axis=1),
            # Count CODING mutations only for quality filtering (aligns with SARS2-mut-fitness)
            # This excludes error-prone noncoding regions from quality criteria while still
            # tracking all mutations for downstream analysis
            n_nt_mutations=lambda x: x["nt_mutations"].map(safe_mutation_count),
            # Count reversions on CODING mutations only
            n_reversions_to_ref=lambda x: x["nt_mutations"].map(
                lambda ms: count_reversions(ms, ref)
            ),
            n_reversions_to_founder=lambda x: x["nt_mutations"].map(
                lambda ms: count_reversions(ms, founder)
            ),
            is_terminal=lambda x: x["leaves_sharing_mutations"] == 1,
            log_size=lambda x: np.log(x["leaves_sharing_mutations"].clip(lower=1)),
        )
        .drop(columns="nt_mutations_coding_and_noncoding")
        .query("n_reversions_to_ref <= @max_reversions_to_ref")
        .query("n_reversions_to_founder <= @max_reversions_to_founder")
        .query("n_nt_mutations <= @max_nt_mutations")
    )

    log(f"Filtered nodes: {len(filtered_nodes)}")

    # Get coding mutations
    coding_mutations = (
        filtered_nodes.drop(columns="noncoding_nt_mutations")
        .explode(["aa_mutations", "nt_mutations", "codon_changes"])
        .assign(
            protein=lambda x: x["aa_mutations"].str.split(":").str[0],
            aa_mutation=lambda x: x["aa_mutations"].str.split(":").str[1],
            synonymous=lambda x: x["aa_mutation"].map(lambda m: m[0] == m[-1]),
        )
        .rename(columns={"nt_mutations": "nt_mutation", "codon_changes": "codon_change"})
        .groupby(["node_id", "nt_mutation", "is_terminal", "log_size"], as_index=False)
        .aggregate(
            protein=pd.NamedAgg("protein", lambda s: ";".join(s)),
            aa_mutation=pd.NamedAgg("aa_mutation", lambda s: ";".join(s)),
            codon_change=pd.NamedAgg("codon_change", lambda s: ";".join(s)),
            synonymous=pd.NamedAgg("synonymous", "all"),
        )
        .assign(noncoding=False)
    )

    log(f"Coding mutations: {len(coding_mutations)}")

    # Get noncoding mutations
    noncoding_mutations = (
        filtered_nodes.drop(columns=["aa_mutations", "nt_mutations", "codon_changes"])
        .explode("noncoding_nt_mutations")
        .rename(columns={"noncoding_nt_mutations": "nt_mutation"})
        .query("nt_mutation.notnull()")
        .assign(
            aa_mutation="noncoding",
            codon_change="noncoding",
            protein="noncoding",
            synonymous=False,
            noncoding=True,
        )[coding_mutations.columns]
    )

    log(f"Noncoding mutations: {len(noncoding_mutations)}")

    # Combine and count
    mutations = pd.concat([coding_mutations, noncoding_mutations], ignore_index=True)

    mutation_counts_tidy = (
        pd.concat(
            [
                mutations.assign(terminal_nodes_only=False),
                mutations.query("is_terminal").assign(terminal_nodes_only=True),
            ]
        )
        .drop(columns="is_terminal")
        .groupby(
            [
                "terminal_nodes_only",
                "protein",
                "aa_mutation",
                "nt_mutation",
                "codon_change",
                "synonymous",
                "noncoding",
            ],
            as_index=False,
            dropna=False,
        )
        .aggregate(
            count=pd.NamedAgg("node_id", "count"),
            mean_log_size=pd.NamedAgg("log_size", "mean"),
        )
        .sort_values("count", ascending=False)
        .assign(
            nt_site=lambda x: x["nt_mutation"].str[1:-1].astype(int),
            reference_nt=lambda x: x["nt_site"].map(lambda r: ref[r - 1]),
            clade_founder_nt=lambda x: x["nt_site"].map(lambda r: founder[r - 1]),
            exclude=lambda x: (
                x["nt_site"].isin(sites_to_exclude)
                | x["nt_mutation"].isin(muts_to_exclude)
                | (x["nt_site"] < site_include_range["start"])
                | (x["nt_site"] > site_include_range["end"])
            ),
        )
    )

    mutation_counts = (
        mutation_counts_tidy.assign(
            terminal_nodes_only=lambda x: x["terminal_nodes_only"].map(
                {False: "count", True: "count_terminal"}
            )
        )
        .pivot_table(
            index=[
                c
                for c in mutation_counts_tidy.columns
                if c not in {"count", "terminal_nodes_only", "mean_log_size"}
            ],
            values="count",
            columns="terminal_nodes_only",
            fill_value=0,
        )
        .assign(
            count=lambda x: x["count"] if "count" in x.columns else 0,
            count_terminal=lambda x: x["count_terminal"]
            if "count_terminal" in x.columns
            else 0,
            count_non_terminal=lambda x: x["count"] - x["count_terminal"],
        )
        .sort_values("count", ascending=False)
        .reset_index()
    )

    # Add log node size
    mutation_counts = mutation_counts.merge(
        mutation_counts_tidy.query("not terminal_nodes_only")[
            ["protein", "aa_mutation", "nt_mutation", "codon_change", "mean_log_size"]
        ],
        on=["protein", "aa_mutation", "nt_mutation", "codon_change"],
        how="outer",
        validate="one_to_one",
    )

mutation_counts.to_csv(output_csv, index=False)
log(f"Counted {len(mutation_counts)} unique mutations for clade {clade}")
log("Done!")
