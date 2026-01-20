"""Computing mutation rates from GLM predictions."""

import pandas as pd
import numpy as np
from . import glm

letters = ["A", "C", "G", "T"]
mut_types = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]


class Rates:
    """Class for computing and storing context-specific mutation rates."""

    def __init__(self, lightswitch_boundaries=None):
        """
        Initialize rate table with all possible mutation conditions.

        Args:
            lightswitch_boundaries: Dict mapping mutation types to genome position
                                    boundaries
        """
        self.lightswitch_boundaries = lightswitch_boundaries or {
            "CT": 13467,
            "AT": 21562,
            "CG": 21562,
            "GC": 21562,
        }

        # Build rate table with all combinations
        conditions = []
        for mut_type in mut_types:
            from_nt = mut_type[0]
            for left_nt in letters:
                for right_nt in letters:
                    motif = left_nt + from_nt + right_nt
                    for unpaired in [0, 1]:
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

        self.rates = pd.DataFrame(conditions)
        self.rates["rate"] = np.nan
        self.rates["predicted_count"] = np.nan
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

        Args:
            df_train: Training DataFrame with curated synonymous mutations
        """
        # Initialize and train GLM
        self.glm = glm.GeneralLinearModel(
            included_factors=["motif", "unpaired", "nt_site_before_boundary"]
        )
        self.glm.train(df_train, self.lightswitch_boundaries)

        # Predict counts for rate table conditions
        self.glm.predict(self.rates)

        # Compute rate from predicted count
        self.rates["rate"] = self.rates["predicted_count"]

    def genome_composition(self, df):
        """
        Calculate genome composition (count of each condition in data).

        Args:
            df: DataFrame with mutation data

        Returns:
            Series of condition counts
        """
        condition_cols = ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
        counts = df.groupby(condition_cols).size().reset_index(name="cond_count")

        self.rates = self.rates.merge(
            counts, on=condition_cols, how="left", suffixes=("", "_new")
        )
        if "cond_count_new" in self.rates.columns:
            self.rates["cond_count"] = (
                self.rates["cond_count_new"].fillna(0).astype(int)
            )
            self.rates.drop(columns=["cond_count_new"], inplace=True)

        return self.rates["cond_count"]

    def predicted_counts_by_clade(self, df):
        """
        Get predicted counts for mutations in a DataFrame.

        Args:
            df: DataFrame with mutations (must have required columns)

        Returns:
            Series of predicted counts
        """
        self.glm.predict(df)
        return df["predicted_count"]

    def residual_variance(self, df, tau):
        """
        Compute residual variance using empirical Bayes shrinkage.

        Args:
            df: DataFrame with mutation data
            tau: Dict of global residual variance by mutation type
        """
        condition_cols = ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]

        # Get empirical residuals per condition
        grouped = df.groupby(condition_cols).apply(
            lambda x: x.actual_count.to_numpy(), include_groups=False
        )
        if len(grouped) == 0:
            return

        idx, emp_counts = zip(*grouped.items())

        idx = pd.MultiIndex.from_tuples(idx, names=condition_cols)
        self.rates = self.rates.set_index(condition_cols)

        emp_residuals = [
            np.mean((np.log(counts + 0.5) - np.log(pred + 0.5)) ** 2)
            for counts, pred in zip(
                emp_counts, self.rates.loc[idx, "predicted_count"].values
            )
        ]
        self.rates.loc[idx, "residual"] = emp_residuals

        # Apply empirical Bayes shrinkage
        av_n_cond = np.mean(self.rates["cond_count"])

        for mut_type in mut_types:
            mask = self.rates.index.get_level_values("mut_type") == mut_type
            n_cond = self.rates.loc[mask, "cond_count"].values
            emp_res = self.rates.loc[mask, "residual"].values

            shrinkage = np.exp(-n_cond / av_n_cond) if av_n_cond > 0 else 1.0
            self.rates.loc[mask, "residual"] = shrinkage * tau.get(
                mut_type, 1.0
            ) + (1 - shrinkage) * emp_res

        self.rates = self.rates.reset_index()
