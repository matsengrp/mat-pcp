#!/usr/bin/env python3
"""Bucket branch mutations into sliding time windows over the full dated tree.

EPIC #115 Phase 2 (#117). Ports the tree-windowing/counting logic from
s2trajectory ``scripts/MATWrapper.py`` (class ``MATWrapper``) into mat-pcp: over
the full filtered tree, group observed branch mutations into overlapping time
windows and, per window, emit each subtree **founder** (with its reconstructed
haplotype) plus per-mutation **actual counts**. No expected counts / no
``delta_fitness_ratio`` (that is Phase 3, #118); no CUSUM / clustering.

Two deliberate deviations from s2trajectory, both locked in the #117 spec:

1. **Windows are anchored at ``end_time`` and tiled backwards** (s2trajectory
   tiles forward from ``start_time``), so the most recent window is always the
   full ``window_size`` and any short/boundary window falls at the *oldest* end.
2. **Branch QC is identical to mat-pcp's clade path** (``count_mutations.py``,
   which inherits the jbloomlab/SARS2-mut-fitness method), not s2trajectory's
   ``qc_filter``: coding mutations only, counted unique (``set``), keep-if
   ``<=`` the three shared config thresholds, plus the unconditional "drop
   branches with >1 nt mutation in one codon" filter. See :func:`qc_filter`.

Design notes:
  * ``import bte`` is **lazy** (inside ``WindowedMAT.__init__``) so the pure
    functions here (``apply_muts``, ``get_codon``, ``qc_filter``,
    ``collapse_multigene``, ...) import cleanly in the base env for the
    ``validate_windowing.py`` checks, which have no bte dependency.
  * Parents are derived from ``.children`` (never ``node.parent`` / ``tree.root``,
    which are unproven in mat-pcp's vendored bte and ``None`` at the root),
    exactly as ``chronumental_qc.py`` does.
"""

from collections import Counter, defaultdict
from datetime import datetime, timedelta
from functools import cache
from warnings import warn

import pandas as pd


# =============================================================================
# Pure helpers (no bte) -- unit-testable in the base env
# =============================================================================

def get_time_dict(chronumental_tsv_file):
    """Map node id -> ``datetime`` from a chronumental ``--dates_out`` TSV.

    Header row, then two columns (node name, predicted date); read positionally
    (``chronumental_qc.py`` does the same), so this does not depend on the column
    labels. Ported from s2trajectory ``MATWrapper.get_time_dict``.
    """
    dates = {}
    with open(chronumental_tsv_file) as fh:
        next(fh)  # skip header
        for line in fh:
            if not line.strip():
                continue
            name, time = line.rstrip("\n").split("\t")[:2]
            dates[name] = datetime.fromisoformat(time)
    return dates


def load_reference_from_fasta(reference_file):
    """Load the single sequence from a one-record FASTA as a string."""
    with open(reference_file) as fh:
        next(fh)  # skip header
        return "".join(map(str.strip, fh))


def apply_muts(sequence, muts):
    """Apply 1-based nucleotide mutations (e.g. ``A25G``) to ``sequence``.

    Warns (does not raise) when a mutation's from-base disagrees with the base
    already present, matching s2trajectory ``MATWrapper.apply_muts``.
    """
    for mut in muts:
        frombase = mut[0]
        tobase = mut[-1]
        site = int(mut[1:-1]) - 1
        if sequence[site] != frombase:
            warn(
                f"apply_muts: from-base {frombase} at site {site + 1} disagrees "
                f"with sequence base {sequence[site]}"
            )
        sequence = sequence[:site] + tobase + sequence[site + 1:]
    return sequence


def get_codon(sequence, start_codon, end_codon, site):
    """Return the 3-nt codon around 1-based ``site`` in ``sequence``.

    The codon offset is inferred from the position at which ``start_codon`` and
    ``end_codon`` differ. Ported verbatim from s2trajectory
    ``MATWrapper.get_codon``.
    """
    offset = [i for i, (a, b) in enumerate(zip(start_codon, end_codon)) if a != b][0]
    return sequence[site - offset - 1: site - offset + 2]


def traverse_mat(root_node, leaf_condition=lambda n: n.is_leaf(), mask=lambda n: False):
    """Preorder traversal from ``root_node`` with custom stop / mask conditions.

    ``leaf_condition(node)`` truthy -> treat as a leaf (stop descending);
    ``mask(child)`` truthy -> skip that child's subtree entirely. Yields
    ``root_node`` even when ``leaf_condition`` is False. Ported from
    s2trajectory. The UShER tree is shallow, so recursion depth is not a concern.
    """
    yield root_node
    if not leaf_condition(root_node):
        for cnode in root_node.children:
            if not mask(cnode):
                yield from traverse_mat(cnode, leaf_condition=leaf_condition, mask=mask)


def iter_leaves(root_node, leaf_condition=lambda n: n.is_leaf(), mask=lambda n: False):
    """Yield only the ``leaf_condition``-satisfying frontier reachable from root.

    Unlike :func:`traverse_mat`, a node with no unmasked children is *not*
    yielded unless ``leaf_condition`` is True on it. Ported from s2trajectory.
    """
    if leaf_condition(root_node):
        yield root_node
    else:
        for cnode in root_node.children:
            if not mask(cnode):
                yield from iter_leaves(cnode, leaf_condition=leaf_condition, mask=mask)


def build_codon_key_map(coding_sites_csv):
    """Build the maps the branch QC needs from ``ref/coding_sites.csv``.

    Returns ``(coding_sites, site_to_codon_keys)`` where:
      * ``coding_sites`` -- set of 1-based sites whose ``gene != "noncoding"``.
      * ``site_to_codon_keys`` -- ``{site: [(gene, codon_site), ...]}`` for each
        coding site. Overlap sites (e.g. ORF1a/ORF1ab) list one ``(gene,
        codon_site)`` pair per gene; the CSV ``;``-joins the per-gene ``gene``
        and ``codon_site`` columns in matching order, so ``zip`` pairs them.

    This is the bte-native equivalent of ``count_mutations.py``'s comma test:
    two distinct nt sites share a codon iff they share a ``(gene, codon_site)``.
    """
    df = pd.read_csv(coding_sites_csv, dtype=str)
    coding_sites = set()
    site_to_codon_keys = {}
    for site_str, codon_site_str, gene_str in zip(
        df["site"], df["codon_site"], df["gene"]
    ):
        if gene_str == "noncoding":
            continue
        site = int(site_str)
        coding_sites.add(site)
        genes = gene_str.split(";")
        codon_sites = codon_site_str.split(";")
        site_to_codon_keys[site] = list(zip(genes, codon_sites))
    return coding_sites, site_to_codon_keys


def count_reversions(mutations, reference_seq):
    """Count unique nt mutations that revert to ``reference_seq``.

    Mirrors ``count_mutations.py:65-76`` (the clade path) so the two analysis
    paths use the same reversion semantics: the reversion test
    (``m[-1] == reference_seq[int(m[1:-1]) - 1]`` over ``set(mutations)``) is
    byte-identical; this copy drops the clade path's ``mutations == [""]``
    sentinel guard and its ``log(...)`` warning (callers here always pass a real
    ``set``, never a ``[""]`` placeholder). (The #117 spec prefers factoring this
    into a shared ``modules/branch_qc.py``; that refactor -- which also touches
    the version-locked clade path -- is deferred, see the PR ``DEFERRED`` section.)
    """
    if not mutations:
        return 0
    try:
        return sum(
            m[-1] == reference_seq[int(m[1:-1]) - 1] for m in set(mutations)
        )
    except (ValueError, IndexError):
        return 0


def qc_filter(
    nt_mutations,
    founder_seq,
    reference_seq,
    coding_sites,
    site_to_codon_keys,
    max_nt_mutations,
    max_reversions_to_ref,
    max_reversions_to_founder,
):
    """Keep-branch predicate, identical to the clade path (``count_mutations.py``).

    Operates on **coding** nt mutations only (noncoding excluded from QC, as in
    ``count_mutations.py:112``), counting **unique** mutations. Keeps the branch
    iff **all four** hold (keep-if ``<=``; criterion 4 unconditional):

      1. ``n_coding_nt_mutations <= max_nt_mutations``
      2. reversions to ``reference_seq`` ``<= max_reversions_to_ref``
      3. reversions to the subtree ``founder_seq`` ``<= max_reversions_to_founder``
      4. no codon carries >1 nt mutation on the branch (SARS2-mut-fitness's
         "ignore any branches with multiple mutations at the same codon";
         ``count_mutations.py:103`` implements this as the comma test).

    ``nt_mutations`` is an iterable of nt-mutation strings (e.g. ``["A25G"]``).
    """
    coding = [m for m in nt_mutations if int(m[1:-1]) in coding_sites]
    unique = set(coding)

    # 1: too many coding mutations
    if len(unique) > max_nt_mutations:
        return False
    # 2: too many reversions to the reference
    if count_reversions(unique, reference_seq) > max_reversions_to_ref:
        return False
    # 3: too many reversions to the subtree founder
    if count_reversions(unique, founder_seq) > max_reversions_to_founder:
        return False
    # 4: >1 nt mutation in one codon -> drop (unconditional)
    codon_to_sites = defaultdict(set)
    for m in unique:
        site = int(m[1:-1])
        for key in site_to_codon_keys.get(site, ()):
            codon_to_sites[key].add(site)
    if any(len(sites) > 1 for sites in codon_to_sites.values()):
        return False

    return True


def collapse_multigene(records, log=None):
    """Collapse per-gene mutation records to mat-pcp ``aggregated.csv`` form.

    ``records`` is an iterable of dicts with keys ``window_start``,
    ``founder_id``, ``protein``, ``nt_site``, ``nt_mutation``, ``aa_mutation``,
    ``codon_change``, ``synonymous``, ``count`` -- one per (gene, nt_mutation).
    A mutation in a gene-overlap region (e.g. ORF1a/ORF1ab) contributes >=2
    records that share the same ``nt_mutation`` and branch ``count``.

    Groups by ``(window_start, founder_id, nt_mutation)`` and, matching
    ``count_mutations.py:144-149``:
      * ``;``-joins ``protein`` / ``aa_mutation`` / ``codon_change`` in
        encounter (row) order -- reproduces ``aggregated.csv``'s ``ORF1a;ORF1ab``
        ordering; do not reorder by config ``genes``.
      * ``count`` = one shared branch count (``max``), never the sum across
        collapsed genes. The per-gene counts are normally equal, but for an
        overlap region a branch can pass the founder-codon guard in one reading
        frame and fail it in the other, so the per-gene subtotals can genuinely
        differ; ``max`` keeps the frame in which the mutation was most observed,
        and a disagreement is logged (via ``log``) rather than silently dropped.
      * ``synonymous`` = ``.all()`` over the collapsed genes (a single bool; the
        ``N``/``ORF9b`` overlap is a distinct reading frame, so disagreement is
        real and ``.all()`` is the intended resolution).

    ``log`` (optional) receives a one-line warning when any collapse group has
    disagreeing per-gene counts. Returns a DataFrame with the fixed 9-column
    ``subtree_counts.csv`` schema.
    """
    columns = [
        "window_start",
        "founder_id",
        "protein",
        "nt_site",
        "nt_mutation",
        "aa_mutation",
        "codon_change",
        "synonymous",
        "count",
    ]
    df = pd.DataFrame(list(records))
    if len(df) == 0:
        return pd.DataFrame(columns=columns)
    group_keys = ["window_start", "founder_id", "nt_mutation"]
    if log is not None:
        # Per-gene counts normally match; a frame-dependent founder-codon guard
        # on overlap genes can make them differ (see docstring). Surface it.
        disagree = df.groupby(group_keys)["count"].nunique()
        n_bad = int((disagree > 1).sum())
        if n_bad:
            log(
                f"collapse_multigene: {n_bad} (window_start, founder_id, "
                f"nt_mutation) group(s) with disagreeing per-gene counts; "
                f"keeping the max"
            )
    collapsed = (
        df.groupby(group_keys, sort=False, as_index=False)
        .agg(
            protein=("protein", lambda s: ";".join(s)),
            nt_site=("nt_site", "first"),
            aa_mutation=("aa_mutation", lambda s: ";".join(s)),
            codon_change=("codon_change", lambda s: ";".join(s)),
            synonymous=("synonymous", "all"),
            count=("count", "max"),
        )
    )
    return collapsed[columns]


# =============================================================================
# WindowedMAT -- the stateful tree wrapper (needs bte)
# =============================================================================

class WindowedMAT:
    """Load a dated MAT once; window over it many times via :meth:`aggregate_data`.

    Ported from s2trajectory ``MATWrapper``. Initialization loads the tree,
    translates mutations, and reconstructs the machinery needed to bucket branch
    mutations into time windows; :meth:`aggregate_data` then runs one window/step
    layout over the pre-loaded state.
    """

    def __init__(
        self,
        tree_path,
        node_dates_path,
        gtf_path,
        fasta_path,
        coding_sites_csv,
        max_nt_mutations,
        max_reversions_to_ref,
        max_reversions_to_founder,
        start_time=None,
        end_time=None,
        log=print,
    ):
        import bte  # lazy: keeps the pure helpers importable without bte

        self.log = log
        self.max_nt_mutations = max_nt_mutations
        self.max_reversions_to_ref = max_reversions_to_ref
        self.max_reversions_to_founder = max_reversions_to_founder

        log(f"Loading tree from {tree_path}")
        self.tree = bte.MATree(tree_path)
        self.nodes = self.tree.depth_first_expansion()
        self.ids_to_nodes = {n.id: n for n in self.nodes}

        # Parents from children (never node.parent / tree.root), as
        # chronumental_qc.py does. parent_of[root] is absent -> .get() is None.
        self.parent_of = {}
        for n in self.nodes:
            for c in n.children:
                self.parent_of[c.id] = n.id
        all_ids = set(self.ids_to_nodes)
        root_ids = all_ids - set(self.parent_of)
        assert len(root_ids) == 1, f"expected exactly one root, found {len(root_ids)}"
        self.root_id = next(iter(root_ids))
        self.root_node = self.ids_to_nodes[self.root_id]

        # Node dates. chronumental dates the full time-tree and Phase 1's
        # chronumental_qc check #1 enforces that every node (root included) is
        # dated, so we assert it here and carry no undated-node contingency.
        self.times = get_time_dict(node_dates_path)
        undated = all_ids - set(self.times)
        assert not undated, (
            f"{len(undated)} tree node(s) have no chronumental date "
            f"(e.g. {sorted(undated)[:5]}); windowing requires all nodes dated"
        )

        if start_time is None:
            self.start_time = self.times[self.root_id]
            log(f"Inferred start_time from root: {self.start_time}")
        else:
            self.start_time = start_time
            log(f"Using configured start_time: {self.start_time}")

        if end_time is None:
            self.end_time = max(self.times.values())
            log(f"Inferred end_time from max node date: {self.end_time}")
        else:
            self.end_time = end_time
            log(f"Using configured end_time: {self.end_time}")

        # Fail fast on a degenerate range rather than silently emitting zero
        # windows (a misconfigured window_start_time/window_end_time, or a
        # single-date tree).
        assert self.start_time < self.end_time, (
            f"start_time ({self.start_time}) must be strictly before end_time "
            f"({self.end_time})"
        )

        self.reference_seq = load_reference_from_fasta(fasta_path)

        log("Translating mutations")
        # {node_id: [bte.AAChange, ...]} for coding mutations; nodes absent from
        # the dict just have no coding parent mutations.
        self.translations = defaultdict(list, self.tree.translate(gtf_path, fasta_path))

        self.coding_sites, self.site_to_codon_keys = build_codon_key_map(
            coding_sites_csv
        )
        log(f"Loaded {len(self.coding_sites)} coding sites")

    @cache
    def get_haplotype(self, node_id):
        # NOTE: functools.cache on a bound method keys on (self, node_id) and so
        # keeps every WindowedMAT instance alive for the process lifetime. Benign
        # here -- each Snakemake script is a short-lived subprocess with a single
        # WindowedMAT -- and inherited from s2trajectory's MATWrapper.
        return self.tree.get_haplotype(node_id)

    def _window_starts(self, window_size, step_size):
        """Window start datetimes, anchored at ``end_time`` and tiled backwards.

        Every window is exactly ``window_size`` wide. ``w_i = end_time -
        window_size - i*step_size``; the newest window ``w_0`` covers
        ``[end_time - window_size, end_time]`` (ends exactly at ``end_time``).
        Emit while the window still overlaps the data range
        (``w_i + window_size > start_time``); the oldest window may extend before
        ``start_time`` (acceptable -- the founder predicate simply admits
        everything from the data start). Returned sorted ascending.
        """
        # Fail fast: a non-positive step (or size) would loop forever below, since
        # ``w`` would never advance past the break condition.
        assert window_size > timedelta(0), f"window_size must be > 0, got {window_size}"
        assert step_size > timedelta(0), f"step_size must be > 0, got {step_size}"
        starts = []
        i = 0
        while True:
            w = self.end_time - window_size - i * step_size
            if w + window_size <= self.start_time:
                break
            starts.append(w)
            i += 1
        starts.sort()
        return starts

    def aggregate_data(self, window_size, step_size):
        """Window the tree once with the given ``window_size`` / ``step_size``.

        Both are ``timedelta``. Returns a list of ``(window_start, subtrees)``
        pairs, one per window (ascending), where each ``subtrees`` entry is the
        :meth:`_summarize_subtree` output:
        ``((founder_id, founder_muts, subtree_size, branches_with_mutations),
        records)`` and each record is
        ``(gene, nt_index, aa_index, aa, nuc, original_codon, alternative_codon,
        is_synonymous, count)``.
        """
        window_starts = self._window_starts(window_size, step_size)
        self.log(
            f"time range {self.start_time} to {self.end_time}; "
            f"window {window_size}, step {step_size}: {len(window_starts)} windows"
        )

        subtree_data = []
        for t in window_starts:
            # Founders of window t: first node inside the window on its lineage.
            # Founder band keys on step_size; subtree extent on window_size (the
            # deliberate step/size asymmetry, s2trajectory 311/340).
            def leaf_func(node, t=t):
                parent_id = self.parent_of.get(node.id)
                return self.times[node.id] >= t and (
                    parent_id is None or self.times[parent_id] < t
                )

            def mask_func(node, t=t):
                return self.times[node.id] > t + step_size

            founders = list(iter_leaves(self.root_node, leaf_func, mask_func))

            subtrees = []
            for founder_node in founders:
                node_list = list(
                    traverse_mat(
                        founder_node,
                        mask=lambda n, t=t: self.times[n.id] > t + window_size,
                    )
                )
                # Keep the subtree iff it has >1 node and >=1 mutation strictly
                # below the founder.
                if len(node_list) > 1 and sum(
                    len(n.mutations) for n in node_list[1:]
                ) > 0:
                    subtrees.append(self._summarize_subtree(node_list))
            subtree_data.append((t, subtrees))
        return subtree_data

    def _summarize_subtree(self, node_list):
        """Aggregate one subtree (founder first) into founder + count records."""
        founder = node_list[0]
        internal_nodes = node_list[1:]
        subtree_size = len(node_list)
        branches_with_mutations = sum(
            1 for n in internal_nodes if len(n.mutations) > 0
        )
        founder_muts = self.get_haplotype(founder.id)
        founder_seq = apply_muts(self.reference_seq, founder_muts)

        records = []
        for node in internal_nodes:
            if not qc_filter(
                node.mutations,
                founder_seq,
                self.reference_seq,
                self.coding_sites,
                self.site_to_codon_keys,
                self.max_nt_mutations,
                self.max_reversions_to_ref,
                self.max_reversions_to_founder,
            ):
                continue
            for aamut in self.translations[node.id]:
                # Recompute the founder codon; skip context-drifted translations
                # whose original_codon no longer matches the founder (s2traj 378).
                founder_codon = get_codon(
                    founder_seq,
                    aamut.original_codon,
                    aamut.alternative_codon,
                    aamut.nt_index,
                )
                if founder_codon == aamut.original_codon:
                    records.append(
                        (
                            aamut.gene,
                            aamut.nt_index,
                            aamut.aa_index,
                            aamut.aa,
                            aamut.nuc,
                            aamut.original_codon,
                            aamut.alternative_codon,
                            aamut.is_synonymous(),
                        )
                    )
        # Dedup + count (Counter preserves insertion / encounter order).
        counted_records = [rec + (count,) for rec, count in Counter(records).items()]
        return (
            (founder.id, founder_muts, subtree_size, branches_with_mutations),
            counted_records,
        )
