#!/usr/bin/env python3
"""Emit one gene's per-clade, per-nucleotide censor table for DASM training.

UShER masking does **not** mark a site ambiguous: the mutations are deleted from
the MAT, so the node inherits its parent's state. A masked site is therefore
*censored, not ambiguous* -- the state is known, the *event* was destroyed. That
distinction is why this is a table of coordinates for the loss rather than an
``N`` written into the PCP sequence: an ``N`` would tokenize the whole codon to
``AA_AMBIG_IDX`` and throw away a state we actually know (#144).

The pipeline's masking output at the time, ``{gene}_codon_sites_to_mask.txt``,
was unusable for this purpose in three independent ways: it never reached the
loss (only the PCP drop filter and the dashboards read it), it was codon-granular
(one suspect nucleotide took its two innocent neighbours with it), and it was
clade-independent by construction -- ``get_problematic_sites.py`` kept only the
reasons that hold in *every* clade, so the entire clade-specific half of
``data/usher_masked_sites.yaml`` was invisible to training. #153 emitted the
clade-independent half at nucleotide granularity and #158 retired the codon file,
so the second and third points are now history; the first is why this table
exists at all.

This table carries exactly the two clade-specific reasons of
``masking.compute_exclusions``' five:

  * ``usher_clade_sites``     -- whole-site masks, blank ``from_nt``/``to_nt``
  * ``usher_clade_mutations`` -- directional masks, both populated

The other three (ProblematicSites VCF, config ``sites_to_exclude``, outside
``site_include_range``) are in ``pcps/{gene}/masked_nt_sites.csv``, at the same
nucleotide granularity since #153, and founder->reference reversions are
deliberately in neither: they
are a selection filter, not an artifact mask. **A training run must apply both
files**, which is why the guards below report the overlap rather than assuming
the two files partition the space.

Coordinates are gene-local, because the consumer has no CDS knowledge -- the PCP
sequence *is* the gene. ``nt_index`` is a 0-based index into that sequence, and
the genomic->index map is one-to-many at a -1 ribosomal frameshift (see
``sequences.genomic_to_pcp_indices``).
"""

import csv
import os
import re
from collections import Counter
from importlib import util as importlib_util

# The by-path module loader is defined once in
# workflow/modules/workflow_loader.py. Because that loader is itself the
# mechanism for importing shared workflow code, it cannot be imported via
# itself -- it is wired up here by absolute path (#111).
_loader_spec = importlib_util.spec_from_file_location(
    "workflow_loader",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "modules",
        "workflow_loader.py",
    ),
)
workflow_loader = importlib_util.module_from_spec(_loader_spec)
_loader_spec.loader.exec_module(workflow_loader)
load_workflow_module = workflow_loader.load_workflow_module

gtf_utils = load_workflow_module("gtf")
masking = load_workflow_module("masking")
sequences = load_workflow_module("sequences")

usher_yaml = snakemake.input.usher_masked_sites
gtf_path = snakemake.input.gtf
ref_fasta = snakemake.input.ref_fasta
output_table = snakemake.output.table
output_roster = snakemake.output.clades
gene = snakemake.wildcards.gene
clades_processed = list(snakemake.params.clades)
founder_fastas = dict(snakemake.params.founder_fastas)
pcp_gene_cds_blocks = snakemake.params.pcp_gene_cds_blocks
log_file = snakemake.log[0] if hasattr(snakemake, "log") and snakemake.log else None

COLUMNS = ["clade", "nt_site", "nt_index", "from_nt", "to_nt"]
MUTATION_RE = re.compile(r"^([ACGT])(\d+)([ACGT])$")


def log(msg):
    """Write to the log file if available, and always to stdout."""
    if log_file:
        with open(log_file, "a") as f:
            f.write(msg + "\n")
    print(msg)


def read_fasta(path):
    """Return a single-record FASTA's sequence, uppercased, newlines stripped."""
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()


def parse_mutation(mutation, clade):
    """Split an ``nt_mutation`` string (``"G26530A"``) into ``(from, site, to)``.

    Raises rather than skipping. A mask the parser cannot read is a mask that
    silently censors nothing, which is the failure mode this whole file exists to
    make impossible -- and it would be invisible in the emitted table, since an
    unparseable entry and an out-of-gene one both contribute zero rows.
    """
    match = MUTATION_RE.match(str(mutation).strip())
    if match is None:
        raise ValueError(
            f"clade {clade!r} has an unparseable UShER mutation mask "
            f"{mutation!r} in {usher_yaml}. Expected nt_mutation format, e.g. "
            f"'G26530A'. Skipping it would silently censor nothing."
        )
    from_nt, site, to_nt = match.group(1), int(match.group(2)), match.group(3)
    if from_nt == to_nt:
        raise ValueError(
            f"clade {clade!r} has UShER mutation mask {mutation!r}, whose parent "
            f"and child bases are identical. That is not a mutation, so it can "
            f"censor nothing; the yaml entry is malformed."
        )
    return from_nt, site, to_nt


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

log(f"Building the censor table for gene: {gene}")

clade_sites, clade_mutations = masking.load_usher_masks(usher_yaml)
log(
    f"Parsed {usher_yaml}: {len(clade_sites)} clades with site masks, "
    f"{sum(1 for c in clade_mutations if clade_mutations[c])} with mutation masks"
)

# resolve_cds_blocks must be called WITH the overrides: ORF1b exists only in
# pcp_gene_cds_blocks, and the bare GTF lookup raises "No CDS blocks found".
cds_blocks = gtf_utils.resolve_cds_blocks(gtf_path, gene, pcp_gene_cds_blocks)
positions = sequences.cds_nt_positions(cds_blocks)
log(f"{gene} CDS blocks {cds_blocks} -> {len(positions)} nt / {len(positions) // 3} codons")

# The roster is part of the deliverable, not an extra. A clade with nothing
# censored contributes zero rows, so in the CSV alone it is indistinguishable
# from a mistyped clade name -- the consumer cannot implement "unknown clade
# raises, empty clade passes" from the rows. Writing the full roster separately
# makes absence-from-roster and absence-from-rows distinguishable.
roster = sorted(clades_processed)
log(f"Clade roster: {len(roster)} clades from the get_clades checkpoint")

# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------

rows = []
sites_outside_gene = 0
muts_outside_gene = 0

for clade in roster:
    for nt_site in sorted(clade_sites.get(clade, ())):
        indices = sequences.genomic_to_pcp_indices(cds_blocks, nt_site)
        if not indices:
            sites_outside_gene += 1
            continue
        # One row per index: a site shared by two blocks at the frameshift is
        # read twice by slice_cds and so occupies two positions of the PCP
        # sequence, in different codons.
        for nt_index in indices:
            rows.append((clade, nt_site, nt_index, "", ""))
    for mutation in sorted(clade_mutations.get(clade, ())):
        from_nt, nt_site, to_nt = parse_mutation(mutation, clade)
        indices = sequences.genomic_to_pcp_indices(cds_blocks, nt_site)
        if not indices:
            muts_outside_gene += 1
            continue
        for nt_index in indices:
            rows.append((clade, nt_site, nt_index, from_nt, to_nt))

# Most yaml masks lie outside any given gene by design, so this is a log line
# rather than an error.
log(
    f"Dropped as outside {gene}: {sites_outside_gene} clade-site masks, "
    f"{muts_outside_gene} clade-mutation masks"
)

# load_usher_masks returns *sets*, which are unordered, so without an explicit
# sort the file is not byte-reproducible across runs -- and data-central
# manifests pin md5/xxh3 per file (CLAUDE.md pattern 9).
rows.sort(key=lambda row: (row[0], row[2], row[4]))

# ---------------------------------------------------------------------------
# Guard 1: the row count reconciles against an independent derivation
# ---------------------------------------------------------------------------
#
# A coordinate error that maps every mask outside the gene writes a well-formed
# empty CSV at exit 0. But emptiness is *legitimate*: config_test_data.yaml's
# max_clades: 3 gives 19A/19B/20A and the yaml's earliest clade is 20B, so the
# CI table is empty for every gene. A "nonzero total" assertion would fail every
# correct CI run. So reconcile against a count computed a different way --
# Counter over the position list, versus the emit loop's per-site block-offset
# search -- and say *why* when the answer is zero.
occurrences = Counter(positions)
expected_rows = 0
for clade in roster:
    for nt_site in clade_sites.get(clade, ()):
        expected_rows += occurrences[nt_site]
    for mutation in clade_mutations.get(clade, ()):
        # Re-parsed rather than carried over from the emit loop on purpose: this
        # count has to be derived independently to be worth anything. Parsing is
        # already known to succeed, since the emit loop ran first.
        expected_rows += occurrences[parse_mutation(mutation, clade)[1]]

if len(rows) != expected_rows:
    raise ValueError(
        f"censor table for {gene!r} has {len(rows)} rows but the UShER masks and "
        f"CDS blocks {cds_blocks} account for {expected_rows}. The emit loop and "
        f"the independent position count disagree, which means a mask was "
        f"dropped or duplicated -- most likely the one-to-many frameshift index."
    )

if not rows:
    clades_in_yaml = set(clade_sites) | set(clade_mutations)
    if not (set(roster) & clades_in_yaml):
        reason = (
            f"none of the {len(roster)} processed clades appears in "
            f"{usher_yaml} (its earliest clade is "
            f"{min(clades_in_yaml) if clades_in_yaml else 'n/a'})"
        )
    else:
        reason = f"no mask of any processed clade falls inside {gene}'s CDS"
    log(f"NOTE: the censor table for {gene} is empty because {reason}. This is "
        f"legitimate; the row count reconciles with the independent derivation.")

# ---------------------------------------------------------------------------
# Guard 2: the clade vocabulary reconciles in both directions
# ---------------------------------------------------------------------------
#
# Reported, never fatal. Both differences have legitimate causes, and the point
# is that a stale or mistyped yaml key becomes visible in the log instead of
# quietly censoring nothing.
in_table = {row[0] for row in rows}
yaml_clades = {c for c in (set(clade_sites) | set(clade_mutations))
               if clade_sites.get(c) or clade_mutations.get(c)}
yaml_not_on_tree = sorted(yaml_clades - set(roster))
roster_unmasked = sorted(set(roster) - in_table)
log(
    f"Clade vocabulary: {len(in_table)} clades contribute rows; "
    f"{len(yaml_not_on_tree)} yaml clades are absent from this tree "
    f"({', '.join(yaml_not_on_tree) or 'none'}); "
    f"{len(roster_unmasked)} processed clades have no mask in {gene}"
)

# ---------------------------------------------------------------------------
# Guard 3: nt_index round-trips, and from_nt is oriented the right way
# ---------------------------------------------------------------------------
#
# The round-trip is stated as a *sequence* identity rather than index
# arithmetic, so it also catches a blocks-vs-sequence mismatch (the wrong gene's
# blocks) that arithmetic alone would pass.
reference = read_fasta(ref_fasta)
gene_reference = sequences.slice_cds(reference, cds_blocks)
if len(gene_reference) != len(positions):
    raise ValueError(
        f"{gene}'s reference CDS slice is {len(gene_reference)} nt but its "
        f"position map has {len(positions)} entries; slice_cds and "
        f"cds_nt_positions disagree about blocks {cds_blocks}."
    )
for clade, nt_site, nt_index, _, _ in rows:
    if gene_reference[nt_index] != reference[nt_site - 1]:
        raise ValueError(
            f"censor row ({clade}, nt_site={nt_site}, nt_index={nt_index}) does "
            f"not round-trip: {gene}'s reference CDS has "
            f"{gene_reference[nt_index]!r} at that index but the genome has "
            f"{reference[nt_site - 1]!r} at that site."
        )

# A directional mask records a *branch* mutation, and matUtils reports each
# branch relative to its parent -- not to the clade founder. So from_nt need not
# equal the founder base: a mask can legitimately sit on a derived background,
# including a reversion back toward the founder. Measured on the v11
# gisaid-and-public tree, of the 1,126 directional masks carried by the 32
# clades that have any (41 of the tree's clades carry a mask of some kind):
# 1,083 have founder == from_nt, 36 have founder == to_nt (a reversion toward
# the founder, exactly the artifact class UShER masking targets), and 7 have
# neither (a twice-derived background, e.g. 25A A22896T where the founder is G).
# Asserting equality -- as #144 originally specified -- would therefore fail on
# real data. What is worth catching is a from_nt/to_nt **swap**, which inverts
# that distribution wholesale, so the ratio is the assertion and the raw counts
# go to the log as the diagnostic.
founder_is_from = founder_is_to = founder_is_neither = 0
founder_cache = {}
for clade, nt_site, nt_index, from_nt, to_nt in rows:
    if not from_nt:
        continue
    if clade not in founder_cache:
        if clade not in founder_fastas:
            raise ValueError(
                f"no {gene} founder slice was passed for clade {clade!r}, which "
                f"the roster names. params.clades and params.founder_fastas are "
                f"built from the same checkpoint call in the Snakefile, so this "
                f"means the two have been decoupled; a bare KeyError here would "
                f"not say that. Known keys: {sorted(founder_fastas)}"
            )
        founder_cache[clade] = read_fasta(founder_fastas[clade])
    founder = founder_cache[clade]
    if len(founder) != len(positions):
        raise ValueError(
            f"{clade}'s {gene} founder slice is {len(founder)} nt but {gene}'s "
            f"CDS is {len(positions)} nt; the founder slice and the censor "
            f"table disagree about this gene's coordinates."
        )
    base = founder[nt_index]
    if base == from_nt:
        founder_is_from += 1
    elif base == to_nt:
        founder_is_to += 1
    else:
        founder_is_neither += 1

directional = founder_is_from + founder_is_to + founder_is_neither
if directional:
    log(
        f"Directional masks vs clade founder: {founder_is_from} have "
        f"founder == from_nt, {founder_is_to} founder == to_nt (reversion "
        f"toward the founder), {founder_is_neither} neither"
    )
    if founder_is_to > founder_is_from:
        raise ValueError(
            f"{founder_is_to} of {directional} directional masks in {gene} have "
            f"the clade founder's base as their *destination* and only "
            f"{founder_is_from} have it as their source. Individually that is "
            f"legitimate (a reversion toward the founder), but a majority is the "
            f"signature of from_nt and to_nt being written the wrong way round, "
            f"which would censor the reverse mutation and leave the real one live."
        )

# ---------------------------------------------------------------------------
# Guard 4: overlap with the global mask is reported, not forbidden
# ---------------------------------------------------------------------------
#
# #144 specified this as a disjointness assertion, on the reasoning that the
# clade-independent mask carries only clade-independent reasons and this table
# only clade-specific ones. The two reason sets are disjoint but the *sites* are
# not: a site can be both clade-masked in the yaml and in the ProblematicSites
# VCF. On the v11 gisaid-and-public tree 32 (clade, nt_site) pairs are both --
# sites 3037, 24378 and 28253 -- so the assertion would fail on real data. The
# consumer's actual requirement is that it applies BOTH files, and for that it
# needs to know the two overlap rather than partition.
#
# Reported as nucleotide positions since #158: through v12 the counterparty was
# codon-granular, so the overlap had to be stated in codons to be comparable at
# all, which understated it by conflating three positions into one.
masked_positions_in_table = sorted({row[2] for row in rows})
log(
    f"{len(masked_positions_in_table)} distinct nt positions of {gene} are "
    f"touched by this table; they may also appear in "
    f"pcps/{gene}/masked_nt_sites.csv, which is expected -- the two files carry "
    f"different reasons, not different sites, and a training run must apply "
    f"both."
)

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

os.makedirs(os.path.dirname(output_table) or ".", exist_ok=True)
with open(output_table, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(COLUMNS)
    writer.writerows(rows)

with open(output_roster, "w") as f:
    f.write("".join(f"{clade}\n" for clade in roster))

log(f"Wrote {len(rows)} rows to {output_table}")
log(f"Wrote {len(roster)} clades to {output_roster}")
log("Done!")
