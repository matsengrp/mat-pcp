"""General Linear Model for mutation rate estimation."""

import numpy as np
import pandas as pd
from scipy import linalg

letters = ["A", "C", "G", "T"]
mut_types = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]


class GeneralLinearModel:
    """GLM for predicting mutation counts based on sequence context."""

    def __init__(self, included_factors, interacting=False, W=None, tau_squared=None):
        """
        Initialize GLM.

        Args:
            included_factors: List of factors to include (e.g., ['motif', 'unpaired'])
            interacting: Whether to include interaction terms
            W: Trained weights (optional, for pre-trained model)
            tau_squared: Residual variance per mutation type (optional)
        """
        self.included_factors = included_factors
        self.interacting = interacting
        self.W = W if W is not None else {}
        self.tau_squared = tau_squared if tau_squared is not None else {}

    def train(self, df_train, lightswitch_boundaries=None):
        """
        Train GLM on mutation counts.

        Args:
            df_train: DataFrame with training data
            lightswitch_boundaries: Dict mapping mutation types to position boundaries
        """
        for mut_type in mut_types:
            nt1, nt2 = mut_type[0], mut_type[1]
            mask = df_train["nt_mutation"].str.match(f"^{nt1}.*{nt2}$").values
            df_mut_type = df_train[mask].copy()

            if len(df_mut_type) == 0:
                self.W[mut_type] = np.zeros((1,))
                self.tau_squared[mut_type] = 1.0
                continue

            # Create design matrix
            X = self.create_data_matrix(df_mut_type, mut_type)

            # Get response variable (log counts with offset)
            y = np.log(df_mut_type["actual_count"].values + 0.5)

            # Fit linear regression
            try:
                self.W[mut_type] = linalg.lstsq(X, y)[0]

                # Compute residual variance
                y_pred = X @ self.W[mut_type]
                residuals = y - y_pred
                self.tau_squared[mut_type] = np.var(residuals)
            except (linalg.LinAlgError, ValueError):
                self.W[mut_type] = np.zeros(X.shape[1])
                self.tau_squared[mut_type] = 1.0

    def create_data_matrix(self, df, mut_type):
        """
        Create design matrix for GLM.

        Args:
            df: DataFrame with mutation data
            mut_type: Mutation type (e.g., 'CT')

        Returns:
            Design matrix X
        """
        n = len(df)

        # Always include intercept
        X = [np.ones(n)]

        if "motif" in self.included_factors:
            # One-hot encode left and right nucleotide context
            context_l = df["motif"].str[0].values
            context_r = df["motif"].str[2].values
            X.extend(self.one_hot_l_r(context_l, context_r))

        if "unpaired" in self.included_factors:
            X.append(df["unpaired"].values.astype(float))

        if "nt_site_before_boundary" in self.included_factors:
            X.append(df["nt_site_before_boundary"].astype(int).values)

        return np.column_stack(X)

    def predict(self, df):
        """
        Add predicted counts to DataFrame.

        Args:
            df: DataFrame with mutation data (modified in place)
        """
        predicted_counts = np.zeros(len(df))
        remaining_variance = np.zeros(len(df))

        for mut_type in mut_types:
            nt1, nt2 = mut_type[0], mut_type[1]
            mask_mut_type = (
                df["nt_mutation"].str.match("^" + nt1 + ".*" + nt2 + "$").values
            )
            df_mut_type = df[mask_mut_type]

            if len(df_mut_type) == 0:
                continue

            X = self.create_data_matrix(df_mut_type.copy(), mut_type)

            if mut_type in self.W and len(self.W[mut_type]) == X.shape[1]:
                predicted_counts[mask_mut_type] = (X @ self.W[mut_type]).flatten()
                remaining_variance[mask_mut_type] = self.tau_squared.get(mut_type, 1.0)
            else:
                predicted_counts[mask_mut_type] = 0
                remaining_variance[mask_mut_type] = 1.0

        df["predicted_count"] = np.exp(predicted_counts) - 0.5
        df["tau_squared"] = remaining_variance

    @staticmethod
    def one_hot_l_r(context_l, context_r):
        """One-hot encode left and right nucleotide context."""
        sigma = {}
        for nt in ["C", "G", "T"]:
            sigma[nt + "_l"] = (context_l == nt).astype(int)
            sigma[nt + "_r"] = (context_r == nt).astype(int)

        return [
            sigma["C_l"],
            sigma["G_l"],
            sigma["T_l"],
            sigma["C_r"],
            sigma["G_r"],
            sigma["T_r"],
        ]
