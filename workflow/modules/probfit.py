"""Module for probabilistic fitness estimation."""

import math
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize


def get_probabilistic_estimate(
    n_obs,
    n_pred,
    tau_squared,
    offset=1.0e-6,
    s_prior=2.5,
    N_f=150,
    N_n_exp=300,
):
    """
    Compute probabilistic fitness estimate using Bayesian inference.

    The model assumes:
    - Observed counts follow a Poisson distribution with rate n_exp * exp(f)
    - n_exp (expected count) follows a log-normal distribution centered on n_pred
    - f (fitness) has a Gaussian prior with mean 0 and std s_prior

    Args:
        n_obs: Observed mutation count
        n_pred: Predicted mutation count from GLM
        tau_squared: Residual variance from GLM
        offset: Small value for numerical stability in log
        s_prior: Standard deviation of Gaussian prior on fitness
        N_f: Number of grid points for fitness integration
        N_n_exp: Number of grid points for expected count integration

    Returns:
        Tuple of (f_hat, left_ci, right_ci, f_mean, f_std)
        - f_hat: Maximum a posteriori fitness estimate
        - left_ci: Left 68% credible interval distance from f_hat
        - right_ci: Right 68% credible interval distance from f_hat
        - f_mean: Posterior mean fitness
        - f_std: Posterior standard deviation
    """
    # Handle edge cases
    if n_pred <= 0:
        return 0.0, 0.0, 0.0, 0.0, 1.0

    n_obs = max(0, int(n_obs))
    tau_squared = max(0.01, tau_squared)

    s2_prior = s_prior * s_prior

    def p_f(f):
        """Prior on fitness."""
        return norm.pdf(f, loc=0.0, scale=s_prior)

    c2 = 1 / math.sqrt(tau_squared * 2 * math.pi)
    c3 = 2 * tau_squared

    def p_n_exp_given_n_pred(n_exp):
        """Log-normal likelihood for expected count."""
        return (
            c2
            * np.exp(-(np.log((n_exp + offset) / n_pred) ** 2) / c3)
            / (n_exp + offset)
        )

    log_n_obs_factorial = math.lgamma(n_obs + 1)

    def p_n_obs_given_f_n_exp(n_exp, f):
        """Poisson likelihood for observed count."""
        log_p = (
            n_obs * (np.log(n_exp + offset) + f)
            - n_exp * np.exp(f)
            - log_n_obs_factorial
        )
        return np.exp(np.clip(log_p, -100, 100))

    # Compute bounds for integration
    mean_log_normal = np.exp(np.log(n_pred) + tau_squared / 2)
    std_log_normal = np.sqrt(
        (np.exp(tau_squared) - 1) * np.exp(2 * np.log(n_pred) + tau_squared)
    )
    upper_bound_log_normal = mean_log_normal + 10 * std_log_normal
    lower_bound_log_normal = max(0, mean_log_normal - 10 * std_log_normal)
    n_obs_plus_1 = 1 + n_obs
    n_obs_plus_1_sqrt = np.sqrt(1 + n_obs)

    def p_n_obs_given_f_n_pred(f):
        """Marginal likelihood integrated over n_exp."""
        e_to_the_minus_f = np.exp(-f)
        mean_poisson = e_to_the_minus_f * n_obs_plus_1
        std_poisson = e_to_the_minus_f * n_obs_plus_1_sqrt
        lower_lim = max(0, max(lower_bound_log_normal, mean_poisson - 5 * std_poisson))
        upper_lim = min(upper_bound_log_normal, mean_poisson + 5 * std_poisson)

        if upper_lim <= lower_lim:
            return offset

        n_exp_values = np.linspace(lower_lim, upper_lim, 100).flatten()
        integrand_values = p_n_exp_given_n_pred(n_exp_values) * p_n_obs_given_f_n_exp(
            n_exp_values, f
        )
        return max(offset, np.trapz(y=integrand_values, x=n_exp_values))

    def p_f_given_n_obs_n_pred(f):
        """Posterior on fitness (unnormalized)."""
        return p_n_obs_given_f_n_pred(f) * p_f(f)

    def get_f_hat():
        """Find MAP estimate of fitness."""
        initial_guess = np.log(n_obs + 0.5) - np.log(n_pred + 0.5)
        result = minimize(
            lambda f: -np.log(p_f_given_n_obs_n_pred(f) + offset),
            x0=initial_guess,
            method="BFGS",
        )
        return result.x[0]

    def get_hdi(x, p, coverage=0.682):
        """Compute highest density interval."""
        dx = x[1] - x[0]
        dp = dx * (p[:-1] + (p[1:] - p[:-1]) / 2)
        sorted_indices = np.argsort(-dp)
        x_sorted = x[sorted_indices]
        dp_sorted = dp[sorted_indices]
        cumulative_prob = np.cumsum(dp_sorted)
        cutoff_idx = np.argmax(cumulative_prob >= coverage)
        hdi_x_values = x_sorted[: cutoff_idx + 1]
        return (min(hdi_x_values), max(hdi_x_values))

    def characterize_posterior(f_max):
        """Characterize posterior distribution."""
        mu_posterior = (
            np.log((n_obs + offset) / n_pred)
            * s2_prior
            / (tau_squared + 1 / (n_obs + offset) + s2_prior)
        )
        sigma_posterior = np.sqrt(
            1 / (1 / (tau_squared + 1 / (n_obs + offset)) + 1 / s2_prior)
        )
        f_lower_lim = mu_posterior - 5 * sigma_posterior
        f_upper_lim = mu_posterior + 5 * sigma_posterior
        f_values = np.linspace(f_lower_lim, f_upper_lim, N_f)

        e_to_the_minus_f_lower = np.exp(-f_lower_lim)
        e_to_the_minus_f_upper = np.exp(-f_upper_lim)
        mean_poisson_lower = e_to_the_minus_f_lower * n_obs_plus_1
        std_poisson_lower = e_to_the_minus_f_lower * n_obs_plus_1_sqrt
        mean_poisson_upper = e_to_the_minus_f_upper * n_obs_plus_1
        std_poisson_upper = e_to_the_minus_f_upper * n_obs_plus_1_sqrt
        n_exp_lower_lim = max(
            0,
            lower_bound_log_normal,
            min(
                mean_poisson_lower - 5 * std_poisson_lower,
                mean_poisson_upper - 5 * std_poisson_upper,
            ),
        )
        n_exp_upper_lim = min(
            upper_bound_log_normal,
            max(
                mean_poisson_lower + 5 * std_poisson_lower,
                mean_poisson_upper + 5 * std_poisson_upper,
            ),
        )

        if n_exp_upper_lim <= n_exp_lower_lim:
            return f_max, 0.1, 0.1, 0.5

        n_exp_values = np.linspace(n_exp_lower_lim, n_exp_upper_lim, N_n_exp)

        n_exp_mesh, f_mesh = np.meshgrid(n_exp_values, f_values, indexing="ij")
        integrand_values = (
            p_n_exp_given_n_pred(n_exp_mesh)
            * p_n_obs_given_f_n_exp(n_exp_mesh, f_mesh)
            * p_f(f_mesh)
        )
        n_exp_integral = np.trapz(integrand_values, x=n_exp_values, axis=0)
        p_n_obs_given_n_pred = np.trapz(n_exp_integral, x=f_values)

        if p_n_obs_given_n_pred <= 0:
            return f_max, 0.1, 0.1, 0.5

        f_mean = (
            np.trapz(f_values * n_exp_integral, x=f_values) / p_n_obs_given_n_pred
        )
        f_std = np.sqrt(
            np.trapz((f_values - f_mean) ** 2 * n_exp_integral, x=f_values)
            / p_n_obs_given_n_pred
        )
        left_interval, right_interval = get_hdi(
            f_values, n_exp_integral / p_n_obs_given_n_pred
        )
        return f_mean, f_max - left_interval, right_interval - f_max, f_std

    f_hat = get_f_hat()
    f_mean, left_ci, right_ci, f_std = characterize_posterior(f_hat)
    return f_hat, left_ci, right_ci, f_mean, f_std


def add_probabilistic_estimates(dataframe, offset=1.0e-6, N_f=150, N_n_exp=300):
    """
    Add probabilistic fitness estimates to a DataFrame.

    Args:
        dataframe: DataFrame with columns actual_count, predicted_count, tau_squared
        offset: Small value for numerical stability
        N_f: Grid points for fitness
        N_n_exp: Grid points for expected counts
    """

    def calculate_estimates(row):
        try:
            return get_probabilistic_estimate(
                row.actual_count,
                row.predicted_count,
                row.tau_squared,
                offset=offset,
                N_f=N_f,
                N_n_exp=N_n_exp,
            )
        except Exception:
            return (0.0, 0.0, 0.0, 0.0, 1.0)

    estimates = dataframe.apply(calculate_estimates, axis=1)
    (
        dataframe["f_max"],
        dataframe["left_conf_int"],
        dataframe["right_conf_int"],
        dataframe["f_mean"],
        dataframe["f_st_dev"],
    ) = zip(*estimates)
