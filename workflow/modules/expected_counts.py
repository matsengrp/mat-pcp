"""Windowed expected counts + ``delta_fitness_ratio`` (EPIC #115 Phase 3, #118).

This module ports s2trajectory's expected-count engine
(``s2trajectory/scripts/ExpectedCalc.py``, class ``PossibleMutations``) into
mat-pcp and drives it with mat-pcp's GLM rates. It is deliberately pure
``pandas``/``numpy`` (no ``bte``): it consumes the Phase-2 sliding-window CSVs
(``windows/founder_muts.csv`` + ``windows/subtree_counts.csv``) together with the
GLM ``master_table.csv``, the reference FASTA, and ``ref/coding_sites.csv``.

Key differences from the s2trajectory original (all deliberate):

* **No ``gene_fix`` hack.** s2trajectory rewrites ORF1a/ORF1b gene labels to work
  around a bug in *its* coding-sites table (its issue #46). mat-pcp's
  ``ref/coding_sites.csv`` is authoritative and is used by *both* the Phase-2
  actual counts and this Phase-3 enumeration, so no relabeling is needed.
* **Rates come from mat-pcp's GLM ``master_table.csv``**, used at their native
  scale (see :func:`prepare_rates_df`) -- no renormalization is needed because the
  per-window ``evol_time`` factor below divides out any constant rate factor.
* **Actual counts join to expected by re-deriving the amino-acid ID from the
  founder enumeration** (merge on ``(site, mut_nt)``), rather than trusting the
  translation in ``subtree_counts.csv``. This guarantees actual and expected sit
  on an identical, founder-context ID system for every window.

Normalization follows **mat-pcp's own clade path** (``rates.py`` ``evol_time``),
NOT s2trajectory's ``rate x subtree_mut_count``: expected counts are scaled
**per window** by ``T_window = (observed synonymous count) / (rate summed over all
possible synonymous mutations that window)``, so total expected synonymous equals
total observed synonymous by construction (see
:func:`compute_fitness_trajectory`). This keeps the windowed fitness on the same
neutral-model footing as the clade-based ``aamut_fitness``. Finally
``delta_fitness_ratio = log((actual + p) / (expected + p))`` with pseudocount
``p`` (default 0.5), matching ``aamutfit.ratio_fitness``.
"""

import re

import numpy as np
import pandas as pd

# Standard genetic code (identical to s2trajectory's CODON_TABLE).
CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #
def load_reference_from_fasta(reference_file):
    """Load the single sequence from a FASTA file and return it as a string."""
    with open(reference_file, "r") as fh:
        next(fh)
        reference_seq = "".join(map(str.strip, fh))
    return reference_seq


def apply_muts(sequence, muts):
    """Return ``sequence`` after applying a ``{site: to_nt}`` mutation dict.

    Sites are 1-based (site 25 is ``sequence[24]``). Ported verbatim from
    s2trajectory so the founder-sequence reconstruction matches.
    """
    if not muts:
        return sequence
    sites = sorted(muts.keys())
    substrings = [sequence[: sites[0] - 1]]
    for start, stop in zip(sites[:-1], sites[1:]):
        substrings.append(muts[start])
        substrings.append(sequence[start : stop - 1])
    substrings.append(muts[sites[-1]])
    substrings.append(sequence[sites[-1] :])
    return "".join(substrings)


def all_poss_muts(site_map):
    """Explode a ``site``/``wt_nt`` map into all three possible alternate nts.

    Returns a tall dataframe with columns ``site``, ``wt_nt``, ``mut_nt``.
    Ported verbatim from s2trajectory.
    """
    ret = site_map.copy()
    alt_nuc = [
        {"A": "C", "T": "C", "C": "T", "G": "C"},
        {"A": "T", "T": "A", "C": "A", "G": "T"},
        {"A": "G", "T": "G", "C": "G", "G": "A"},
    ]
    re_pattern = re.compile("A|T|C|G")
    for j in range(3):
        fn = lambda x, j=j: alt_nuc[j][x.group(0)]
        ret[j + 1] = ret.wt_nt.str.replace(re_pattern, fn, regex=True)
    ret_tall = ret.melt(id_vars=["site", "wt_nt"], value_name="mut_nt").drop(
        "variable", axis=1
    )
    return ret_tall


def wide_codon_sites_to_tall(coding_sites):
    """Expand semicolon-joined overlap-gene coding-site rows into one row per gene.

    Input: a dataframe with columns ``site``, ``codon_position``, ``codon_site``,
    ``gene`` where overlapping genes are ``;``-joined within a single row (the
    ``ref/coding_sites.csv`` format). Output: one row per (site, gene) with columns
    ``site``, ``gene``, ``codon_position``, ``codon_sites``.
    """
    records = {"site": [], "gene": [], "codon_position": [], "codon_sites": []}
    for _, data in coding_sites.iterrows():
        genes = str(data.gene).split(";")
        codon_poss = str(data.codon_position).split(";")
        codon_sites = str(data.codon_site).split(";")
        for gene, codon_pos, codon_site in zip(genes, codon_poss, codon_sites):
            records["site"].append(data.site)
            records["gene"].append(gene)
            records["codon_position"].append(codon_pos)
            records["codon_sites"].append(codon_site)
    return pd.DataFrame(records)


def parse_mutation_string(mut):
    """Parse a nt-mutation string like ``A25G`` into ``(from_nt, site, to_nt)``."""
    return mut[0], int(mut[1:-1]), mut[-1]


def founder_muts_to_dict(founder_nt_mutations, bad_sites):
    """Convert a space-joined founder haplotype into a ``{site: to_nt}`` dict.

    Sites in ``bad_sites`` are dropped (matching s2trajectory's driver, which
    filters founder mutations by the bad-site list before reconstruction).
    """
    if not isinstance(founder_nt_mutations, str) or not founder_nt_mutations.strip():
        return {}
    out = {}
    for mut in founder_nt_mutations.split():
        _, site, to_nt = parse_mutation_string(mut)
        if site not in bad_sites:
            out[site] = to_nt
    return out


# --------------------------------------------------------------------------- #
# Input adapters (mat-pcp -> s2trajectory-style tables)
# --------------------------------------------------------------------------- #
def build_site_unpaired_map(rna_structure_file, genome_len):
    """Build a per-site ``unpaired`` (bool) table from the Lan 2022 structure file.

    Replicates ``annotate_counts.py``'s ``assign_ss_pred`` so the ``unpaired``
    feature used here matches the one the GLM was trained on: a site is
    ``unpaired`` iff its pairing partner (5th whitespace field) is 0; sites absent
    from the file default to paired (``unpaired = False``), matching the GLM's
    ``"nd" -> 0`` behavior.

    Returns a dataframe with columns ``site`` (1..genome_len) and ``unpaired``
    (bool).
    """
    with open(rna_structure_file) as f:
        lines = [line.rstrip().split() for line in f]
    paired = np.array([[int(x[0]), int(x[4])] for x in lines[1:]])
    paired_dict = dict(zip(paired[:, 0], paired[:, 1]))
    sites = np.arange(1, genome_len + 1)
    unpaired = np.array([paired_dict.get(int(s), None) == 0 for s in sites], dtype=bool)
    return pd.DataFrame({"site": sites, "unpaired": unpaired})


def prepare_rates_df(master_table_csv):
    """Load mat-pcp's GLM per-context rates from ``glm/master_table.csv``.

    Keeps the per-context columns the founder-context merge needs
    (``mut_type``, ``motif``, ``unpaired``, ``nt_site_boundary``,
    ``nt_site_before_boundary``) plus mat-pcp's ``rate`` column **as-is**.

    The rate's absolute scale does not matter here: expected counts are
    normalized per window by an ``evol_time``-style factor ``T`` (see
    :func:`compute_fitness_trajectory`), exactly as mat-pcp's clade path does
    (``rates.py``). ``T`` divides out any constant factor on ``rate`` (including
    mat-pcp's ``n_sites`` factor), so no fraction-scale renormalization is needed.
    ``expected_count`` / ``residual`` are dropped; ``unpaired`` and
    ``nt_site_before_boundary`` are cast to bool to match the merge dtypes.
    """
    master = pd.read_csv(master_table_csv)
    rates_df = master[
        [
            "mut_type",
            "motif",
            "unpaired",
            "nt_site_boundary",
            "nt_site_before_boundary",
            "rate",
        ]
    ].copy()
    rates_df["unpaired"] = rates_df["unpaired"].astype(bool)
    rates_df["nt_site_before_boundary"] = rates_df["nt_site_before_boundary"].astype(
        bool
    )
    return rates_df


def parse_id_metadata(ids):
    """Recover per-mutation metadata from the amino-acid ``ID`` strings.

    ``ID`` is ``<gene>-<wt_aa><aa_site><mut_aa>`` (gene names carry no ``-`` and
    amino acids are single characters, ``*`` for stop), so gene/site/aa are fully
    recoverable without carrying the (large) per-founder tables around. Returns a
    dataframe with columns ``ID``, ``gene``, ``aa_site``, ``aa_mutation``,
    ``wt_aa``, ``mut_aa``, ``synonymous``.
    """
    rows = []
    for the_id in ids:
        gene, aa = the_id.split("-", 1)
        wt_aa, mut_aa, aa_site = aa[0], aa[-1], int(aa[1:-1])
        rows.append((the_id, gene, aa_site, aa, wt_aa, mut_aa, wt_aa == mut_aa))
    return pd.DataFrame(
        rows,
        columns=[
            "ID", "gene", "aa_site", "aa_mutation", "wt_aa", "mut_aa", "synonymous",
        ],
    )


# --------------------------------------------------------------------------- #
# Expected-count engine (port of PossibleMutations)
# --------------------------------------------------------------------------- #
class ExpectedCounts:
    """Compute per-subtree expected counts from founder sequence context.

    Parameters mirror s2trajectory's ``PossibleMutations`` but take already-loaded
    objects instead of file paths, and drop the ``gene_fix`` hack.
    """

    def __init__(self, ref_seq, coding_sites_csv, bad_sites, basepair_df, rates_df):
        self.ref_seq = ref_seq
        self.bad_site_list = set(bad_sites)
        self.basepair_df = basepair_df  # columns: site, unpaired(bool)
        self.rates_df = rates_df.copy()
        # Per-mut_type site boundary, then drop it (only its per-mut_type value is
        # needed to recompute nt_site_before_boundary in founder context).
        self.boundary_df = self.rates_df.groupby("mut_type").nt_site_boundary.first()
        self.rates_df.drop(columns=["nt_site_boundary"], inplace=True)

        coding_sites = pd.read_csv(coding_sites_csv)
        if len(coding_sites) != len(self.ref_seq):
            raise ValueError(
                "reference sequence and coding sites must be the same length, got "
                f"{len(self.ref_seq)} and {len(coding_sites)}"
            )
        if not np.all(coding_sites.site == np.arange(1, len(self.ref_seq) + 1)):
            raise ValueError("coding sites must be in order of the reference sequence")
        self.coding_sites = coding_sites
        self.tall_coding_sites = wide_codon_sites_to_tall(coding_sites)

        self.reference_site_map = pd.DataFrame(
            {"site": np.arange(1, len(self.ref_seq) + 1), "wt_nt": list(self.ref_seq)}
        )
        self.reference_site_map = self.reference_site_map.query(
            "~site.isin(@self.bad_site_list)"
        )
        self._keep_sites = set(self.reference_site_map.site)
        # The reference enumeration + rates is computed once; each founder only
        # recomputes the codons its mutations touch (s2trajectory's optimization).
        # Founder tables are NOT cached: subtree founders are almost always unique
        # per (founder, window), so a cache would give no reuse while retaining a
        # ~90k-row table per founder (multi-GB on large trees).
        self._reference_table = self._build_rate_table(
            self.reference_site_map, self.ref_seq
        )

    def _possible_mutations_df(self, site_map, ref_seq):
        """All possible coding mutations for ``site_map`` under ``ref_seq``.

        Returns a dataframe with ``site``, ``wt_nt``, ``mut_nt``, ``gene``,
        ``codon_position``, ``codon_sites``, ``wt_codon``, ``mut_codon``, ``wt_aa``,
        ``mut_aa``, ``is_synonymous``, ``nuc``, ``aa``, ``ID``.
        """
        annotated = self.tall_coding_sites.query(
            "site.isin(@site_map.site.to_list())"
        )
        all_muts = all_poss_muts(site_map)
        all_muts = all_muts.merge(annotated, on="site", how="outer")

        # Drop noncoding rows (they produce NaN codon info), then cast to int.
        all_muts.drop(
            index=all_muts.query("codon_position == 'noncoding'").index, inplace=True
        )
        all_muts.codon_position = all_muts.codon_position.astype(int)

        all_muts["wt_codon"] = (all_muts.site - all_muts.codon_position).apply(
            lambda x: ref_seq[x : x + 3]
        )
        mut_codons = []
        for j in range(1, 4):
            rows = all_muts.query("codon_position == @j")
            slc = rows.wt_codon.str.slice
            mut_codons.append(slc(stop=j - 1) + rows.mut_nt + slc(start=j))
        all_muts["mut_codon"] = pd.concat(mut_codons)
        all_muts["wt_aa"] = all_muts.wt_codon.map(CODON_TABLE)
        all_muts["mut_aa"] = all_muts.mut_codon.map(CODON_TABLE)
        all_muts["is_synonymous"] = all_muts.wt_aa == all_muts.mut_aa
        all_muts["nuc"] = all_muts.wt_nt + all_muts.site.astype(str) + all_muts.mut_nt
        all_muts["aa"] = all_muts.wt_aa + all_muts.codon_sites.astype(str) + all_muts.mut_aa
        all_muts["ID"] = all_muts.gene + "-" + all_muts.aa
        return all_muts

    def _add_context_features(self, the_df, ref_seq):
        """Add ``mut_type``, ``motif``, ``unpaired``, ``nt_site_before_boundary``.

        Left-joins (not inner) so a site/mut_type missing from ``basepair_df`` /
        ``boundary_df`` surfaces as NaN and fails fast, rather than silently
        dropping the mutation.
        """
        motif_at_site = lambda s: ref_seq[s - 2 : s + 1]
        the_df = the_df.merge(self.basepair_df, on="site", how="left", validate="m:1")
        the_df["motif"] = the_df.site.apply(motif_at_site)
        the_df["mut_type"] = the_df.wt_nt + the_df.mut_nt
        the_df = the_df.merge(
            self.boundary_df, on="mut_type", how="left", validate="m:1"
        )
        missing = the_df["unpaired"].isna() | the_df["nt_site_boundary"].isna()
        if missing.any():
            raise ValueError(
                f"{int(missing.sum())} mutation rows lack an unpaired/boundary "
                "annotation (site or mut_type not covered); this should be complete."
            )
        the_df["nt_site_before_boundary"] = the_df.site < the_df.nt_site_boundary
        return the_df

    def _build_rate_table(self, site_map, ref_seq):
        """Enumerate + context-annotate + rate-merge one ``site_map``/``ref_seq``.

        Returns a table with columns ``ID``, ``gene``, ``aa_site``, ``wt_aa``,
        ``mut_aa``, ``is_synonymous``, ``site``, ``mut_nt``, ``aa_mutation``,
        ``rate`` (mat-pcp's GLM per-context rate for each mutation's context under
        ``ref_seq``). Raises if any row fails to match a GLM rate.
        """
        poss = self._possible_mutations_df(site_map, ref_seq)
        poss = self._add_context_features(poss, ref_seq)
        poss = poss.merge(
            self.rates_df,
            on=["mut_type", "motif", "unpaired", "nt_site_before_boundary"],
            how="left",
            validate="m:1",
        )
        n_missing = int(poss["rate"].isna().sum())
        if n_missing:
            # The GLM rate grid is complete (every mut_type x motif x unpaired x
            # boundary combination), so a missing rate signals a real bug -- a
            # context/dtype mismatch in the merge. Fail fast rather than silently
            # contaminating expected counts with zero rates.
            bad = poss.loc[
                poss["rate"].isna(),
                ["site", "mut_type", "motif", "unpaired", "nt_site_before_boundary"],
            ].head(5)
            raise ValueError(
                f"{n_missing} possible-mutation rows had no matching GLM rate; the "
                f"master-table context grid should be complete. First few:\n{bad}"
            )
        table = poss[
            [
                "ID", "gene", "codon_sites", "wt_aa", "mut_aa", "is_synonymous",
                "site", "mut_nt", "aa", "rate",
            ]
        ].rename(columns={"codon_sites": "aa_site", "aa": "aa_mutation"})
        # Drop ORF1a in favor of ORF1ab (matching aamutfit.py's convention): in
        # mat-pcp's coding_sites the entire ORF1a region is labeled "ORF1a;ORF1ab",
        # so every such mutation is enumerated under both genes with an identical
        # codon/frame -- keeping only ORF1ab removes that redundant duplicate and
        # aligns gene labels with the clade-based aamut_fitness output.
        table = table[table["gene"] != "ORF1a"]
        return table

    def _touched_sites(self, founder_muts):
        """Sites whose GLM rate can change under this founder (bad sites removed).

        A founder mutation shifts the rate of two distinct site sets, and the
        recompute must cover **both** or a subset of rates goes stale:

        1. **Translation** -- a mutation rewrites its whole codon, so every site
           in that codon needs re-translating (``wt_aa``/``is_synonymous``).
        2. **Motif context** -- the GLM rate keys on the 3-nt motif
           ``ref_seq[s-2:s+1]`` *centered* on ``s``, so a mutation at site ``X``
           changes the motif of its immediate neighbors ``X-1`` and ``X+1`` even
           when they lie in a *different* codon. This ±1 halo is keyed on every
           founder mutation (coding or not -- a noncoding mutation still shifts a
           coding neighbor's motif).

        Recompute the union of both. (The halo was previously omitted, silently
        staling the rate of coding sites adjacent to a founder mutation across a
        codon boundary -- roughly one neighbor per codon-boundary mutation.)
        """
        mut_sites = list(founder_muts.keys())
        codon_data = self.tall_coding_sites.query(
            "site.isin(@mut_sites) & codon_position != 'noncoding'"
        )
        starts = codon_data.site - codon_data.codon_position.astype(int) + 1
        sites = {start + j for start in starts for j in range(3)}
        sites.update(x + d for x in mut_sites for d in (-1, 0, 1))
        return sorted(sites & self._keep_sites)

    def founder_rate_table(self, founder_muts):
        """Return the per-mutation GLM rate table for one founder subtree.

        ``founder_muts`` is a ``{site: to_nt}`` dict (bad sites already removed).
        Built by splicing recomputed rows for the codons the founder mutations
        touch into the shared reference enumeration -- so only a handful of codons
        are re-translated per founder. Columns: ``ID``, ``gene``, ``aa_site``,
        ``wt_aa``, ``mut_aa``, ``is_synonymous``, ``site``, ``mut_nt``,
        ``aa_mutation``, ``rate``.
        """
        if not founder_muts:
            return self._reference_table

        sites = self._touched_sites(founder_muts)
        if not sites:
            return self._reference_table

        founder_seq = apply_muts(self.ref_seq, founder_muts)
        wt_nts = [founder_muts.get(s, self.ref_seq[s - 1]) for s in sites]
        founder_site_map = pd.DataFrame({"site": sites, "wt_nt": wt_nts})
        recomputed = self._build_rate_table(founder_site_map, founder_seq)

        touched = set(sites)
        table = pd.concat(
            [
                self._reference_table[~self._reference_table.site.isin(touched)],
                recomputed,
            ],
            ignore_index=True,
        )
        return table


# --------------------------------------------------------------------------- #
# Driver: windowed actual + expected -> delta_fitness_ratio trajectory
# --------------------------------------------------------------------------- #
def compute_fitness_trajectory(
    founder_muts_df,
    subtree_counts_df,
    ec,
    bad_sites,
    pseudocount=0.5,
    min_subtree_size=0,
    log=print,
):
    """Build the mutation x time-window ``delta_fitness_ratio`` trajectory.

    Expected counts use mat-pcp's ``evol_time`` normalization (``rates.py``), the
    same neutral-model application as the clade-based path, computed **per window**:

        T_window = (observed synonymous count in the window)
                   / (sum of GLM rate over all possible synonymous mutations across
                      the window's subtree founders)
        expected(ID) = T_window * sum_over_founders rate(context of ID)

    ``T`` divides out the rate's absolute scale, so mat-pcp's raw ``rate`` is used
    directly, and ``sum(expected_syn) == sum(actual_syn)`` holds per window by
    construction (the calibration that makes synonymous ``delta_fitness_ratio``
    center at 0). Then ``delta_fitness_ratio = log((actual + p) / (expected + p))``.

    Args:
        founder_muts_df: ``windows/founder_muts.csv`` as a dataframe.
        subtree_counts_df: ``windows/subtree_counts.csv`` as a dataframe.
        ec: an initialized :class:`ExpectedCounts`.
        bad_sites: set of globally excluded sites.
        pseudocount: fitness pseudocount ``p`` in ``log((a + p) / (e + p))``.
        min_subtree_size: drop founders whose ``subtree_size`` is ``<=`` this
            (s2trajectory's ``min_founder_children`` gate; default 0 keeps all).
        log: logging callable.

    Returns a tall dataframe with columns ``window_start``, ``ID``, ``gene``,
    ``aa_site``, ``aa_mutation``, ``wt_aa``, ``mut_aa``, ``synonymous``,
    ``actual_count``, ``expected_count``, ``delta_fitness_ratio``.
    """
    founders = founder_muts_df.copy()
    n_before = len(founders)
    founders = founders[founders["subtree_size"] > min_subtree_size]
    if len(founders) < n_before:
        log(
            f"Dropped {n_before - len(founders)} founder-windows with "
            f"subtree_size <= {min_subtree_size}"
        )
    founders["founder_muts"] = founders["founder_nt_mutations"].apply(
        lambda s: founder_muts_to_dict(s, bad_sites)
    )

    unmatched_obs = 0
    total_obs = 0
    skipped_windows = 0
    out_frames = []

    for window_start, w_founders in founders.groupby("window_start", sort=True):
        rate_by_id = None   # sum over founders of GLM rate, grouped by amino-acid ID
        act_by_id = None    # sum over founders of observed counts, grouped by ID
        syn_ids = set()     # amino-acid IDs that are synonymous (wt_aa == mut_aa)
        w_counts = subtree_counts_df[subtree_counts_df["window_start"] == window_start]

        for _, frow in w_founders.iterrows():
            founder_id = frow["founder_id"]
            fmuts = frow["founder_muts"]
            frt = ec.founder_rate_table(fmuts)
            syn_ids.update(frt.loc[frt["is_synonymous"], "ID"].unique())

            contrib = frt.groupby("ID")["rate"].sum()
            rate_by_id = (
                contrib if rate_by_id is None
                else rate_by_id.add(contrib, fill_value=0.0)
            )

            # Actual per ID: attach each observed (site, mut_nt) count to the
            # founder-context amino-acid ID(s) at that site (broadcasts to each
            # overlap gene). The observed wt_nt is intentionally ignored -- every
            # observed mutation is interpreted in the founder's sequence context, so
            # actual and expected share one ID system.
            fc = w_counts[w_counts["founder_id"] == founder_id]
            if len(fc):
                obs = fc[["nt_site", "nt_mutation", "count"]].copy()
                obs["mut_nt"] = obs["nt_mutation"].str[-1]
                obs = obs.rename(columns={"nt_site": "site"})
                total_obs += len(obs)
                merged = obs.merge(
                    frt[["site", "mut_nt", "ID"]], on=["site", "mut_nt"], how="left"
                )
                unmatched_obs += int(merged["ID"].isna().sum())
                merged = merged.dropna(subset=["ID"])
                if len(merged):
                    act_contrib = merged.groupby("ID")["count"].sum()
                    act_by_id = (
                        act_contrib if act_by_id is None
                        else act_by_id.add(act_contrib, fill_value=0.0)
                    )

        if rate_by_id is None:
            continue

        window_df = rate_by_id.rename("rate_by_id").to_frame()
        window_df["actual_count"] = (
            act_by_id if act_by_id is not None else 0.0
        )
        window_df["actual_count"] = window_df["actual_count"].fillna(0.0)
        # Index-alignment safety: every observed ID must be a possible-mutation ID.
        if act_by_id is not None:
            assert set(act_by_id.index) <= set(window_df.index), (
                "observed-count IDs are not a subset of enumerated IDs"
            )

        # evol_time-style per-window normalization (rates.py evol_time): T calibrates
        # so that total expected synonymous == total observed synonymous this window.
        syn_mask = window_df.index.isin(syn_ids)
        syn_rate_total = float(window_df.loc[syn_mask, "rate_by_id"].sum())
        syn_actual_total = float(window_df.loc[syn_mask, "actual_count"].sum())
        if syn_rate_total <= 0 or syn_actual_total == 0:
            skipped_windows += 1
            log(
                f"Skipping window {window_start}: no usable synonymous calibration "
                f"(syn_rate_total={syn_rate_total:.4g}, "
                f"syn_actual_total={syn_actual_total:.0f})"
            )
            continue
        t_window = syn_actual_total / syn_rate_total

        window_df["expected_count"] = t_window * window_df["rate_by_id"]
        window_df["window_start"] = window_start
        out_frames.append(
            window_df.drop(columns=["rate_by_id"]).reset_index()
        )

    if total_obs:
        log(
            f"Observed count rows: {total_obs}; unmatched to founder enumeration: "
            f"{unmatched_obs} ({100.0 * unmatched_obs / total_obs:.2f}%)"
        )
    if skipped_windows:
        log(f"Skipped {skipped_windows} window(s) lacking synonymous calibration data")

    if not out_frames:
        cols = [
            "window_start", "ID", "gene", "aa_site", "aa_mutation", "wt_aa",
            "mut_aa", "synonymous", "actual_count", "expected_count",
            "delta_fitness_ratio",
        ]
        return pd.DataFrame(columns=cols)

    result = pd.concat(out_frames, ignore_index=True)
    meta = parse_id_metadata(result["ID"].unique())
    result = result.merge(meta, on="ID", how="left")
    result["delta_fitness_ratio"] = np.log(
        (result["actual_count"] + pseudocount)
        / (result["expected_count"] + pseudocount)
    )
    result = result[
        [
            "window_start", "ID", "gene", "aa_site", "aa_mutation", "wt_aa",
            "mut_aa", "synonymous", "actual_count", "expected_count",
            "delta_fitness_ratio",
        ]
    ].sort_values(["window_start", "gene", "aa_site", "mut_aa"]).reset_index(drop=True)
    return result
