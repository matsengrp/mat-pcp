"""Computing mutation rates from GLM predictions.

Synchronized with SARS2-mut-fitness-v2/modules/rates.py to ensure consistent
fitness estimates between pipelines.
"""

import os
import sys
import pandas as pd
import numpy as np
from importlib import util as importlib_util


def _load_glm_module():
    """Load the glm module from the same directory."""
    # Get the directory containing this file
    module_dir = os.path.dirname(os.path.abspath(__file__))
    glm_path = os.path.join(module_dir, "glm.py")

    if "glm" in sys.modules:
        return sys.modules["glm"]

    spec = importlib_util.spec_from_file_location("glm", glm_path)
    glm_module = importlib_util.module_from_spec(spec)
    sys.modules["glm"] = glm_module
    spec.loader.exec_module(glm_module)
    return glm_module


glm = _load_glm_module()

letters = ["A", "C", "G", "T"]
mut_types = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]

# Default light-switch boundaries for SARS-CoV-2 mutation types
# These match the boundaries used in SARS2-mut-fitness-v2
DEFAULT_LIGHTSWITCH_BOUNDARIES = {
    "CT": 13467,  # End of ORF1ab
    "AT": 21562,  # End of S gene
    "CG": 21562,
    "GC": 21562,
    # The other eight of the twelve mut_types -- AC, AG, CA, GA, GT, TA, TC, TG --
    # are intentionally omitted: no boundary feature. (This comment previously
    # listed only six, omitting CA and GT; the count is pinned by
    # tests/test_lightswitch_boundary.py.)
}


class Rates:
    """Class for computing and storing context-specific mutation rates."""

    def __init__(self, lightswitch_boundaries=None):
        """
        Initialize rate table with all possible mutation conditions.

        Args:
            lightswitch_boundaries: Dict mapping mutation types to genome position
                boundaries. If None, uses DEFAULT_LIGHTSWITCH_BOUNDARIES.
        """
        self.lightswitch_boundaries = (
            lightswitch_boundaries
            if lightswitch_boundaries is not None
            else DEFAULT_LIGHTSWITCH_BOUNDARIES.copy()
        )

        # Build rate table with all combinations
        # Only include nt_site_before_boundary=True rows for mutation types with boundaries
        conditions = []
        for mut_type in mut_types:
            from_nt = mut_type[0]
            for left_nt in letters:
                for right_nt in letters:
                    motif = left_nt + from_nt + right_nt
                    for unpaired in [0, 1]:
                        # For mutation types WITH boundaries: create both True and False rows
                        # For mutation types WITHOUT boundaries: only create False rows
                        if mut_type in self.lightswitch_boundaries:
                            for before_boundary in [False, True]:
                                conditions.append(
                                    {
                                        "mut_type": mut_type,
                                        "motif": motif,
                                        "unpaired": unpaired,
                                        "nt_site_before_boundary": before_boundary,
                                        "condition": f"{mut_type}_{motif}_"
                                        f"{unpaired}_{before_boundary}",
                                    }
                                )
                        else:
                            # No boundary for this mutation type - only False row
                            conditions.append(
                                {
                                    "mut_type": mut_type,
                                    "motif": motif,
                                    "unpaired": unpaired,
                                    "nt_site_before_boundary": False,
                                    "condition": f"{mut_type}_{motif}_"
                                    f"{unpaired}_False",
                                }
                            )

        self.rates = pd.DataFrame(conditions)
        self.rates["rate"] = np.nan
        self.rates["expected_count"] = np.nan
        self.rates["cond_count"] = 0
        self.rates["residual"] = np.nan

        # Create nt_mutation column for GLM predict
        self.rates["nt_mutation"] = (
            self.rates["mut_type"].str[0] + "1" + self.rates["mut_type"].str[1]
        )

        self.glm = None

    def populate_rates(self, df_train):
        """
        Train GLM and populate rate table with predictions.

        Matches SARS2-mut-fitness-v2 approach:
        1. Train GLM on training data
        2. Predict rates for all conditions
        3. Rescale by total mutations / number of sites

        Args:
            df_train: Training DataFrame with curated synonymous mutations
        """
        # Initialize and train GLM (uses L2 regularization matching v2)
        self.glm = glm.GeneralLinearModel(
            included_factors=["global_context", "rna_structure", "local_context"],
            lightswitch_boundaries=self.lightswitch_boundaries
        )
        self.glm.train(df_train, self.lightswitch_boundaries)

        # Add temporary nt_site column for GLM prediction (matches SARS2-mut-fitness-v2)
        # For nt_site_before_boundary=False, use max site; for True, use 0
        s_max = df_train["nt_site"].max() if "nt_site" in df_train.columns else 30000
        self.rates["nt_site"] = self.rates["nt_site_before_boundary"].apply(
            lambda x: 0 if x else s_max
        )

        # Predict counts for rate table conditions
        self.glm.predict(self.rates, self.lightswitch_boundaries)

        # Drop temporary nt_site column
        self.rates.drop(columns=["nt_site"], inplace=True)

        # Store raw predicted count as rate
        self.rates["rate"] = self.rates["expected_count"]

        # Rescale rates by total number of mutations and number of sites
        # (matches v2 populate_rates normalization)
        tot_mut = df_train["actual_count"].sum()
        n_sites = df_train.shape[0]

        if tot_mut > 0:
            self.rates["rate"] = self.rates["rate"] * n_sites / tot_mut

        # Validate: warn about significantly negative rates
        # Note: Sparse data can cause negative rates in GLM predictions
        # Original issue #55 had TA/TG rates of -0.003 to -0.004
        NEGATIVE_RATE_TOLERANCE = -1e-6
        negative_rates = self.rates[self.rates["rate"] < NEGATIVE_RATE_TOLERANCE]
        if len(negative_rates) > 0:
            # Log warning but don't fail (sparse test data can have unstable predictions)
            import sys
            negative_summary = negative_rates[
                ["mut_type", "motif", "unpaired", "nt_site_before_boundary", "rate"]
            ].head(10)
            print(f"WARNING: {len(negative_rates)} conditions have negative rates (threshold: {NEGATIVE_RATE_TOLERANCE})",
                  file=sys.stderr)
            print(f"First 10 negative rates:\n{negative_summary.to_string()}", file=sys.stderr)
            print("Note: This can occur with sparse training data", file=sys.stderr)

        # Add condition tuple column for lookup
        cond_cols = ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
        self.rates["condition"] = self.rates[cond_cols].apply(tuple, axis=1)

    def genome_composition(self, df):
        """
        Calculate genome composition (count of each condition in data).

        Matches SARS2-mut-fitness-v2 approach - returns array of counts
        aligned with self.rates rows.

        Args:
            df: DataFrame with mutation data

        Returns:
            np.array of condition counts aligned with self.rates
        """
        gen_comp = np.array(
            self.rates.apply(
                lambda x: sum(
                    (x.mut_type == df.mut_type)
                    & (x.motif == df.motif)
                    & (x.unpaired == df.unpaired)
                    & (x.nt_site_before_boundary == df.nt_site_before_boundary)
                ),
                axis=1,
            )
        )
        return gen_comp

    def evol_time(self, df_syn):
        """
        Compute evolutionary time normalization factor T.

        This is the key clade-specific normalization that ensures predicted
        counts are scaled appropriately for each clade's mutation count.

        T = N_mutations / (genome_composition @ rates)

        Args:
            df_syn: DataFrame with synonymous mutations for a clade

        Returns:
            float: Normalization factor T
        """
        gen_comp = self.genome_composition(df_syn)
        N_mut = df_syn["actual_count"].sum()

        denominator = gen_comp @ np.array(self.rates["rate"])
        if denominator > 0:
            T = N_mut / denominator
        else:
            T = 1.0

        return T

    def expected_counts_for_clade(self, df_syn):
        """
        Update rates with clade-specific predicted counts and genome composition.

        Matches SARS2-mut-fitness-v2 predicted_counts() method.

        Args:
            df_syn: DataFrame with synonymous mutations for a clade
        """
        T = self.evol_time(df_syn)
        gen_comp = self.genome_composition(df_syn)

        assert self.rates.shape[0] == len(gen_comp)

        self.rates["expected_count"] = T * self.rates["rate"]
        self.rates["cond_count"] = gen_comp

    def expected_counts_by_clade(self, df):
        """
        Get predicted counts for mutations in a DataFrame using rate lookup.

        Matches SARS2-mut-fitness-v2 approach - looks up predicted counts
        from the rate table by condition tuple.

        Args:
            df: DataFrame with mutations (must have required columns)

        Returns:
            np.array of predicted counts
        """
        # Look-up dictionary: condition -> predicted count
        count_dict = self.rates.set_index("condition")["expected_count"].to_dict()

        # List of row-wise conditions
        cond_list = np.column_stack(
            [
                df["mut_type"],
                df["motif"],
                df["unpaired"],
                df["nt_site_before_boundary"],
            ]
        )

        # Vector of predicted counts via lookup
        pred_count = np.array(
            list(map(lambda x: count_dict.get(tuple(x), 0.0), cond_list))
        )

        return pred_count

    def residual_variance(self, df, tau):
        """
        Compute residual variance using empirical Bayes shrinkage.

        Matches SARS2-mut-fitness-v2 implementation exactly:
        1. Initialize residuals to 0 (not NaN)
        2. Compute empirical residuals for conditions with data
        3. Apply shrinkage: when cond_count=0, shrinkage=1, so residual=tau

        Args:
            df: DataFrame with mutation data
            tau: Dict of global residual variance by mutation type
        """
        # Initialize residuals to 0 (matches v2)
        self.rates["residual"] = np.zeros(self.rates.shape[0], float)

        # Get indices of conditions with data
        idx = self.rates[self.rates["cond_count"] != 0].index

        # Group by condition to get actual counts
        emp_counts = (
            df.groupby(["mut_type", "motif", "unpaired", "nt_site_before_boundary"])
            .apply(lambda x: x.actual_count.to_numpy(), include_groups=False)
            .to_list()
        )

        if len(idx) != len(emp_counts):
            return

        # Compute empirical residuals
        LOG_COUNT_OFFSET = glm.LOG_COUNT_OFFSET
        self.rates.loc[idx, "residual"] = list(
            map(
                lambda x, y: np.mean((np.log(x + LOG_COUNT_OFFSET) - np.log(y + LOG_COUNT_OFFSET)) ** 2),
                emp_counts,
                self.rates.loc[idx, "expected_count"].to_numpy(),
            )
        )

        # Apply empirical Bayes shrinkage (matches v2)
        av_n_cond = np.mean(self.rates["cond_count"])

        self.rates["residual"] = (
            self.rates.groupby("mut_type")
            .apply(
                lambda x: np.exp(-x.cond_count / av_n_cond) * tau[x.name]
                + (1 - np.exp(-x.cond_count / av_n_cond)) * x.residual,
                include_groups=False,
            )
            .explode()
            .values
        )

    def residual_by_clade(self, df):
        """
        Get residuals for mutations in a DataFrame using rate lookup.

        Matches SARS2-mut-fitness-v2 approach.

        Args:
            df: DataFrame with mutations (must have required columns)

        Returns:
            np.array of residuals (tau_squared values)
        """
        # Look-up dictionary: condition -> residual
        res_dict = self.rates.set_index("condition")["residual"].to_dict()

        # List of row-wise conditions
        cond_list = np.column_stack(
            [
                df["mut_type"],
                df["motif"],
                df["unpaired"],
                df["nt_site_before_boundary"],
            ]
        )

        # Vector of residuals via lookup
        res = np.array(
            list(map(lambda x: res_dict.get(tuple(x), 0.1), cond_list))
        )

        return res


#: Era key for the single pooled fit. :func:`add_expected_count` is
#: :func:`add_expected_count_by_era` with one era under this key, so the pooled
#: path is not a separate implementation that could drift from the per-era one.
#: This is what makes ISSUE-137 criterion 5 ("single-era degeneracy is exact")
#: true by construction rather than by measurement.
POOLED = "pooled"


def _fit_rates(train_df, lightswitch_boundaries):
    """Fit one :class:`Rates` on ``train_df`` and return it with its site set.

    The whole fit for one training population: GLM, per-mutation-type residual
    variance ``tau``, genome composition, and empirical-Bayes shrunk residuals.
    Factored out of :func:`add_expected_count` so a per-era fit is the *same*
    code as the pooled fit rather than a copy of it.

    ``train_df`` is mutated (an ``expected_count`` column is added), matching the
    pre-#137 behaviour of :func:`add_expected_count`.

    Returns:
        ``(rate, sites)`` -- the fitted ``Rates`` and the ``nt_site`` values it
        was trained on. The site set is returned rather than recomputed by the
        caller because ``T`` must be normalized over *this* fit's training sites:
        with two eras there are two site sets, and using the wrong one would
        rescale a clade's expected counts by a plausible-looking factor.
    """
    rate = Rates(lightswitch_boundaries)

    sites = train_df["nt_site"].unique()

    # Train GLM on training data
    rate.populate_rates(train_df)

    # Compute initial predicted counts for training data
    train_df["expected_count"] = rate.expected_counts_by_clade(train_df)

    # Compute tau (residual variance) per mutation type from training data
    LOG_COUNT_OFFSET = glm.LOG_COUNT_OFFSET
    tau = train_df.groupby("mut_type").apply(
        lambda x: np.mean(
            (np.log(x.actual_count + LOG_COUNT_OFFSET) - np.log(x.expected_count + LOG_COUNT_OFFSET)) ** 2
        ),
        include_groups=False,
    )

    # Update genome composition from training data
    rate.rates["cond_count"] = rate.genome_composition(train_df)

    # Compute residual variance with empirical Bayes shrinkage
    rate.residual_variance(train_df, tau.to_dict())

    return rate, sites


def add_expected_count_by_era(
    train_by_era,
    count_df,
    clades,
    era_of_clade,
    lightswitch_boundaries=None,
    log=None,
):
    """Add expected counts to ``count_df``, fitting one model per era (#137).

    ``SARS2-mut-fitness-v2`` fits separate neutral models for pre-Omicron and
    Omicron clades; mat-pcp lost that covariate in the port. This restores it.

    Args:
        train_by_era: ``{era: training DataFrame}``. One fit per key. Keys are
            opaque here -- this function never inspects them, so the era
            *definition* stays in ``modules/era.py`` and this stays testable
            against synthetic era labels.
        count_df: counts for all clades; ``expected_count`` and ``tau_squared``
            are written in place.
        clades: clade names to process.
        era_of_clade: callable mapping a clade to its key in ``train_by_era``.
            Pass ``era.era_for_clade``, which raises on an unmapped clade --
            **never** a function with a default era. Upstream dispatches on the
            positional slice ``clades[0:11]``, correct only while the counts file
            happens to arrive in Nextstrain clade order; dispatching on the label
            is the point of this function.
        lightswitch_boundaries: passed to each :class:`Rates`.
        log: optional ``callable(str)`` for progress. Which eras were fit, with
            how many training rows each, is reported here because a fit on the
            wrong partition produces a plausible rate table rather than an error.

    Raises:
        ValueError: if a clade's era has no fitted model. This is deliberately
            *not* caught alongside the empty-era skip: an absent era is decided
            at DAG-construction time by requesting only the non-empty eras'
            tables (``Snakefile: eras_present``), so reaching this point means a
            clade was routed to an era that was never fit -- a wiring bug, not a
            legitimate skip. One ``try/except`` around both turns a broken fit
            into a silently absent table.
    """
    if log is None:
        def log(_msg):
            pass

    if not train_by_era:
        raise ValueError(
            "train_by_era is empty: no era has training data, so no rates can "
            "be fitted. An era with no clades is skipped at DAG-construction "
            "time; an empty mapping here means every era was skipped."
        )

    fits = {}
    for era in sorted(train_by_era):
        train_df = train_by_era[era]
        log(f"Fitting era {era!r} on {len(train_df)} curated training rows")
        fits[era] = _fit_rates(train_df, lightswitch_boundaries)
        rate, sites = fits[era]
        log(
            f"  era {era!r}: {len(sites)} training sites, rate range "
            f"{rate.rates['rate'].min():.6g} - {rate.rates['rate'].max():.6g}"
        )

    # Process each clade against its own era's fit.
    for c in clades:
        era = era_of_clade(c)
        if era not in fits:
            raise ValueError(
                f"clade {c!r} maps to era {era!r}, which has no fitted rate "
                f"model (fitted eras: {sorted(fits)}). A clade must be routed to "
                f"an era that was actually fit; an era present in the counts "
                f"cannot be empty, so this is a wiring error rather than the "
                f"legitimate empty-era skip."
            )
        rate, sites = fits[era]

        count_clade = count_df.loc[count_df["clade"] == c, :].copy()
        count_clade_syn = count_clade.loc[
            count_clade["mut_class"] == "synonymous"
        ].copy()
        # Retain only sites present in training data
        count_clade_syn = count_clade_syn.loc[count_clade_syn["nt_site"].isin(sites)]

        # Compute clade-specific normalization
        rate.expected_counts_for_clade(count_clade_syn)

        # Get predicted counts and residuals for all mutations in clade
        pred_count_clade = rate.expected_counts_by_clade(count_clade)
        residual_clade = rate.residual_by_clade(count_clade)

        # Update count_df with predicted counts and tau_squared
        count_df.loc[count_clade.index, "expected_count"] = pred_count_clade
        count_df.loc[count_clade.index, "tau_squared"] = residual_clade


def add_expected_count(train_df, count_df, clades, lightswitch_boundaries=None):
    """
    Add expected counts to count_df for all clades, from one pooled fit.

    This function matches the SARS2-mut-fitness-v2 add_predicted_count function:
    1. Train GLM on training data
    2. For each clade, compute clade-specific normalization factor T
    3. Apply normalized expectations to all mutations in the clade

    Since #137 this is :func:`add_expected_count_by_era` with a single era, so
    the pooled and per-era paths cannot diverge: the pooled result *is* the
    degenerate per-era result, which is ISSUE-137 criterion 5 by construction.

    Args:
        train_df: Training DataFrame with curated synonymous mutations
        count_df: DataFrame with mutation counts for all clades
        clades: List of clade names to process
        lightswitch_boundaries: Dict mapping mutation types to position boundaries
    """
    add_expected_count_by_era(
        train_by_era={POOLED: train_df},
        count_df=count_df,
        clades=clades,
        era_of_clade=lambda _clade: POOLED,
        lightswitch_boundaries=lightswitch_boundaries,
    )
