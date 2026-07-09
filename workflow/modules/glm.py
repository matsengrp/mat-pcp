"""General Linear Model for mutation rate estimation.

Synchronized with SARS2-mut-fitness-v2/modules/glm.py to ensure consistent
fitness estimates between pipelines.
"""

import numpy as np
import pandas as pd

letters = ["A", "C", "G", "T"]
mut_types = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]

# L2 regularization strength (matches SARS2-mut-fitness-v2)
REGULARIZATION_STRENGTH = 0.1

# Offset for log transformation to handle zero counts in GLM
# Using log(count + 0.5) transformation to:
# 1. Handle zero counts (log(0) is undefined)
# 2. Stabilize variance for Poisson-distributed count data
# 3. Enable linear regression on log-transformed scale
LOG_COUNT_OFFSET = 0.5


class GeneralLinearModel:
    """GLM for predicting mutation counts based on sequence context."""

    def __init__(self, included_factors, interacting=False, W=None, tau_squared=None, lightswitch_boundaries=None):
        """
        Initialize GLM.

        Args:
            included_factors: List of factors to include (e.g., ['global_context', 'rna_structure', 'local_context'])
            interacting: Whether to include interaction terms
            W: Trained weights (optional, for pre-trained model)
            tau_squared: Residual variance per mutation type (optional)
            lightswitch_boundaries: Dict mapping mutation types to genome position boundaries (optional)
        """
        self.included_factors = included_factors
        self.interacting = interacting
        self.W = W if W is not None else {}
        self.tau_squared = tau_squared if tau_squared is not None else {}
        self.lightswitch_boundaries = lightswitch_boundaries if lightswitch_boundaries is not None else {}

    def train(self, df_train, lightswitch_boundaries=None):
        """
        Train GLM on mutation counts using L2-regularized linear regression.

        Uses the same approach as SARS2-mut-fitness-v2:
        - L2 regularization with strength 0.1
        - Intercept term is NOT regularized
        - Tau squared computed as mean squared error on log scale

        Args:
            df_train: DataFrame with training data
            lightswitch_boundaries: Dict mapping mutation types to position boundaries
        """
        # Store boundaries for use during prediction
        if lightswitch_boundaries is not None:
            self.lightswitch_boundaries = lightswitch_boundaries
        # Check which mutation types are present in training data
        present_mut_types = (
            df_train["nt_mutation"].apply(lambda x: x[0] + x[-1]).unique()
        )

        for mut_type in present_mut_types:
            nt1, nt2 = mut_type[0], mut_type[1]
            mask = df_train["nt_mutation"].str.match(f"^{nt1}.*{nt2}$").values
            df_mut_type = df_train[mask].copy()

            if len(df_mut_type) == 0:
                self.W[mut_type] = np.zeros((1,))
                self.tau_squared[mut_type] = 1.0
                continue

            # Create design matrix
            X = self.create_data_matrix(df_mut_type, mut_type, self.lightswitch_boundaries)

            # Get response variable (log counts with offset) - shape (n, 1)
            log_counts = np.log(df_mut_type["actual_count"].values + LOG_COUNT_OFFSET).reshape(-1, 1)

            # Fit L2-regularized linear regression (matches v2 exactly)
            try:
                regularization_matrix = REGULARIZATION_STRENGTH * np.identity(
                    X.shape[1]
                )
                regularization_matrix[0, 0] = 0  # Don't regularize intercept
                w = (
                    np.linalg.inv(X.T @ X + regularization_matrix) @ X.T @ log_counts
                )
                self.W[mut_type] = w

                # Compute residual variance as mean squared error (matches v2)
                self.tau_squared[mut_type] = np.mean((log_counts - X @ w) ** 2)
            except np.linalg.LinAlgError:
                self.W[mut_type] = np.zeros((X.shape[1], 1))
                self.tau_squared[mut_type] = 1.0

    def _compute_global_context_feature(self, df, mut_type, lightswitch_boundaries):
        """
        Compute global context (boundary) feature for a mutation type.

        Args:
            df: DataFrame with mutation data
            mut_type: Mutation type (e.g., 'CT')
            lightswitch_boundaries: Dict mapping mutation types to position boundaries

        Returns:
            np.ndarray or None: Binary array indicating positions after boundary,
                               or None if mutation type has no configured boundary

        Raises:
            ValueError: If mut_type has a configured boundary but 'nt_site' column is missing
        """
        if mut_type not in lightswitch_boundaries:
            return None

        if "nt_site" not in df.columns:
            raise ValueError(
                f"Cannot compute global_context for {mut_type}: 'nt_site' column missing. "
                f"Boundary configured at position {lightswitch_boundaries[mut_type]}."
            )

        boundary = lightswitch_boundaries[mut_type]
        return (df["nt_site"] > boundary).values.astype(int)

    def _compute_local_context_features(self, df):
        """
        Extract and one-hot encode left/right nucleotide context from motif.

        Args:
            df: DataFrame with 'motif' column containing 3-character strings

        Returns:
            list of np.ndarray: One-hot encoded context features
        """
        left_context = df["motif"].str[0].values
        right_context = df["motif"].str[2].values
        return self.one_hot_l_r(left_context, right_context)

    def create_data_matrix(self, df, mut_type, lightswitch_boundaries=None):
        """
        Create design matrix for GLM.

        Maps DataFrame columns to GLM factor names:
        - DataFrame 'motif' column → GLM 'local_context' factor (one-hot encoded left/right nucleotides)
        - DataFrame 'unpaired' column → GLM 'rna_structure' factor (binary RNA structure indicator)
        - DataFrame 'nt_site' column → GLM 'global_context' factor (binary: site > boundary position)

        The GLM uses 'included_factors' with these internal factor names, while
        DataFrames and output CSVs continue using the original column names for compatibility.

        Args:
            df: DataFrame with mutation data. Must include:
                - 'motif': 3-character string (e.g., 'ACA')
                - 'unpaired': 0/1 for RNA structure (if rna_structure factor included)
                - 'nt_site': Genome position (required if mut_type has configured boundary)
            mut_type: Mutation type (e.g., 'CT', 'TA')
            lightswitch_boundaries: Dict mapping mutation types to position boundaries.
                Only mutation types present in this dict will have the global_context
                feature included. If None or empty, no boundary features are added.

        Returns:
            Design matrix X (np.ndarray)

        Raises:
            ValueError: If mut_type has a configured boundary but 'nt_site' column is missing
        """
        if lightswitch_boundaries is None:
            lightswitch_boundaries = {}

        factor_cols = {}

        # Get global context - ONLY for mutation types with configured boundaries
        global_context = self._compute_global_context_feature(
            df, mut_type, lightswitch_boundaries
        )
        factor_cols["global_context"] = [global_context] if global_context is not None else None

        # Get RNA structure (from 'unpaired' column)
        factor_cols["rna_structure"] = [df["unpaired"].values] if "unpaired" in df.columns else None

        # Get left and right context (from 'motif' column)
        factor_cols["local_context"] = self._compute_local_context_features(df)

        # Offset term (intercept)
        base = np.full(len(df), 1)

        # Add columns depending on which factors are included
        columns = [base]
        for factor in self.included_factors:
            if factor_cols[factor] is not None:  # KEY: Skip factors set to None
                columns += factor_cols[factor]

        # Compile data matrix
        X = np.column_stack(columns)
        return X

    def predict(self, df, lightswitch_boundaries=None):
        """
        Add predicted counts to DataFrame.

        Args:
            df: DataFrame with mutation data (modified in place)
            lightswitch_boundaries: Dict mapping mutation types to position boundaries.
                If None, uses the boundaries stored during training.
        """
        # Use stored boundaries if not provided
        if lightswitch_boundaries is None:
            lightswitch_boundaries = self.lightswitch_boundaries

        # Check which mutation types are present in the data
        present_mut_types = df["nt_mutation"].apply(lambda x: x[0] + x[-1]).unique()

        expected_counts = np.zeros(len(df))
        remaining_variance = np.zeros(len(df))

        for mut_type in present_mut_types:
            nt1, nt2 = mut_type[0], mut_type[1]
            mask_mut_type = (
                df["nt_mutation"].str.match("^" + nt1 + ".*" + nt2 + "$").values
            )
            df_mut_type = df[mask_mut_type]

            if len(df_mut_type) == 0:
                continue

            X = self.create_data_matrix(df_mut_type.copy(), mut_type, lightswitch_boundaries)

            # Fail-fast: validate weights exist
            if mut_type not in self.W:
                raise ValueError(
                    f"No trained weights found for mutation type {mut_type}. "
                    f"Model must be trained before making predictions. "
                    f"Available mutation types: {list(self.W.keys())}"
                )

            w = self.W[mut_type]
            # Handle both column vector (n, 1) and 1D array (n,) formats
            if w.ndim == 2:
                w = w.flatten()

            # Fail-fast: validate weight dimensions match design matrix
            if len(w) != X.shape[1]:
                raise ValueError(
                    f"Weight dimension mismatch for {mut_type}: "
                    f"weights have {len(w)} parameters but design matrix has {X.shape[1]} columns. "
                    f"This indicates the model was trained with different features or boundaries."
                )

            expected_counts[mask_mut_type] = (X @ w).flatten()
            remaining_variance[mask_mut_type] = self.tau_squared.get(
                mut_type, 1.0
            )

        df["expected_count"] = np.exp(expected_counts) - LOG_COUNT_OFFSET
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
