"""
LRD v6 — Core Metrics Module

This module implements:
- DFA-2 (Detrended Fluctuation Analysis with quadratic detrending)
- Block bootstrap confidence intervals
- Phase-randomized surrogate tests

Design Principles:
- Self-contained (NumPy/SciPy only for core functionality)
- Conservative methodology following Peng et al. (1994), Kantelhardt et al. (2002)
- Explicit uncertainty quantification

Author: Igor Chechelnitsky
License: CC BY 4.0
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Callable

import numpy as np
from numpy.fft import rfft, irfft
from scipy import stats


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class DFAResult:
    """Results from Detrended Fluctuation Analysis."""
    alpha: float              # Scaling exponent
    intercept: float          # Log-log intercept
    r2: float                 # Fit quality (R²)
    scales: np.ndarray        # Window sizes used
    fluct: np.ndarray         # Fluctuation values F(s)
    fit_range: Tuple[int, int]  # Indices used for fitting
    n_points: int = 0         # Original series length
    
    def __post_init__(self):
        if self.n_points == 0:
            self.n_points = len(self.fluct)


@dataclass
class BootstrapResult:
    """Results from block bootstrap confidence interval estimation."""
    point_estimate: float
    ci_lower: float
    ci_upper: float
    confidence_level: float
    n_bootstrap: int
    block_size: int
    bootstrap_distribution: np.ndarray = field(default_factory=lambda: np.array([]))


@dataclass
class SurrogateResult:
    """Results from phase-randomized surrogate test."""
    observed_value: float
    surrogate_mean: float
    surrogate_std: float
    p_value_two_sided: float
    n_surrogates: int
    surrogate_distribution: np.ndarray = field(default_factory=lambda: np.array([]))


@dataclass
class LRDAnalysisResult:
    """Complete LRD analysis results for a time series."""
    domain: str
    series_name: str
    n: int
    
    # DFA results
    alpha: float
    H: float
    r2: float
    
    # Bootstrap CI
    ci_lower: float
    ci_upper: float
    
    # Surrogate test
    surrogate_mean: float
    surrogate_std: float
    surrogate_p: float
    
    # Quality flags
    is_significant: bool = False  # p < 0.05 for surrogate test
    is_stable: bool = False       # CI doesn't include 0.5 for volatility, etc.
    
    def to_dict(self) -> Dict:
        return {
            "domain": self.domain,
            "series_name": self.series_name,
            "n": self.n,
            "alpha": round(self.alpha, 4),
            "H": round(self.H, 4),
            "r2": round(self.r2, 4),
            "ci_lower": round(self.ci_lower, 4),
            "ci_upper": round(self.ci_upper, 4),
            "surrogate_mean": round(self.surrogate_mean, 4),
            "surrogate_std": round(self.surrogate_std, 4),
            "surrogate_p": round(self.surrogate_p, 4),
            "is_significant": self.is_significant,
            "is_stable": self.is_stable
        }


# =============================================================================
# Utility Functions
# =============================================================================

def _clean_series(x: np.ndarray) -> np.ndarray:
    """Remove NaN and Inf values from series."""
    x = np.asarray(x, dtype=np.float64)
    mask = np.isfinite(x)
    if not mask.any():
        raise ValueError("Series contains no finite values")
    return x[mask]


def _validate_series(x: np.ndarray, min_length: int = 128) -> np.ndarray:
    """Validate and clean time series for analysis."""
    x = _clean_series(x)
    if len(x) < min_length:
        raise ValueError(f"Series too short: {len(x)} < {min_length}")
    return x


# =============================================================================
# DFA-2 Implementation
# =============================================================================

def dfa(
    x: np.ndarray,
    order: int = 2,
    min_scale: int = 16,
    max_scale: Optional[int] = None,
    n_scales: int = 30,
    fit_range: Optional[Tuple[int, int]] = None
) -> DFAResult:
    """
    Detrended Fluctuation Analysis (DFA).
    
    Parameters
    ----------
    x : np.ndarray
        Input time series (will be cleaned of NaN/Inf)
    order : int
        Polynomial detrending order (2 = DFA-2, recommended)
    min_scale : int
        Minimum window size (default: 16)
    max_scale : int or None
        Maximum window size (default: len(x) // 4)
    n_scales : int
        Number of logarithmically spaced scales
    fit_range : tuple or None
        Indices (start, end) into scale array for fitting.
        If None, uses central 60% of scales.
    
    Returns
    -------
    DFAResult
        Contains alpha (scaling exponent), r2, scales, fluctuations
    
    References
    ----------
    Peng C-K et al. (1994). Quantification of scaling exponents and 
    crossover phenomena in nonstationary heartbeat time series. Chaos.
    
    Kantelhardt JW et al. (2002). Multifractal detrended fluctuation 
    analysis of nonstationary time series. Physica A.
    """
    x = _validate_series(x, min_length=min_scale * 4)
    n = len(x)
    
    if max_scale is None:
        max_scale = max(min_scale + 1, n // 4)
    
    # Step 1: Build integrated profile
    # Y(k) = Σᵢ₌₁ᵏ (xᵢ - ⟨x⟩)
    y = np.cumsum(x - np.mean(x))
    
    # Step 2: Generate logarithmically spaced scales
    scales = np.unique(np.floor(
        np.logspace(np.log10(min_scale), np.log10(max_scale), n_scales)
    ).astype(int))
    scales = scales[scales >= min_scale]
    
    if len(scales) < 6:
        raise ValueError(
            f"Not enough scales: {len(scales)}. "
            "Increase series length or adjust scale parameters."
        )
    
    # Step 3: Calculate fluctuations for each scale
    fluct = np.zeros(len(scales), dtype=np.float64)
    
    for i, s in enumerate(scales):
        n_segments = n // s
        if n_segments < 2:
            fluct[i] = np.nan
            continue
        
        F2_values = []
        
        # Forward segments
        for seg in range(n_segments):
            idx_start = seg * s
            idx_end = idx_start + s
            segment_y = y[idx_start:idx_end]
            t = np.arange(s)
            
            # Polynomial fit (detrending)
            coeffs = np.polyfit(t, segment_y, order)
            trend = np.polyval(coeffs, t)
            residual = segment_y - trend
            F2_values.append(np.mean(residual ** 2))
        
        # Backward segments (reduces boundary effects)
        for seg in range(n_segments):
            idx_end = n - seg * s
            idx_start = idx_end - s
            if idx_start < 0:
                continue
            segment_y = y[idx_start:idx_end]
            t = np.arange(s)
            
            coeffs = np.polyfit(t, segment_y, order)
            trend = np.polyval(coeffs, t)
            residual = segment_y - trend
            F2_values.append(np.mean(residual ** 2))
        
        fluct[i] = np.sqrt(np.mean(F2_values))
    
    # Step 4: Log-log linear fit
    valid = np.isfinite(fluct) & (fluct > 0)
    scales_valid = scales[valid]
    fluct_valid = fluct[valid]
    
    if len(scales_valid) < 4:
        raise ValueError("Not enough valid fluctuation points for fitting")
    
    # Determine fit range
    if fit_range is None:
        # Default: central 60% of scales
        i0 = int(0.2 * len(scales_valid))
        i1 = int(0.8 * len(scales_valid))
    else:
        i0, i1 = fit_range
    
    i0 = max(0, min(i0, len(scales_valid) - 3))
    i1 = max(i0 + 2, min(i1, len(scales_valid)))
    
    log_scales = np.log(scales_valid[i0:i1])
    log_fluct = np.log(fluct_valid[i0:i1])
    
    slope, intercept, r_value, _, _ = stats.linregress(log_scales, log_fluct)
    
    return DFAResult(
        alpha=float(slope),
        intercept=float(intercept),
        r2=float(r_value ** 2),
        scales=scales_valid,
        fluct=fluct_valid,
        fit_range=(i0, i1),
        n_points=n
    )


# =============================================================================
# Block Bootstrap
# =============================================================================

def block_bootstrap_ci(
    x: np.ndarray,
    estimator: Callable[[np.ndarray], float],
    n_bootstrap: int = 500,
    block_size: Optional[int] = None,
    confidence: float = 0.95,
    seed: Optional[int] = None
) -> BootstrapResult:
    """
    Block bootstrap confidence interval for time series estimator.
    
    Uses non-overlapping block resampling to preserve temporal
    dependence structure within blocks.
    
    Parameters
    ----------
    x : np.ndarray
        Input time series
    estimator : callable
        Function that takes array and returns scalar estimate
    n_bootstrap : int
        Number of bootstrap replications
    block_size : int or None
        Block length (default: sqrt(n))
    confidence : float
        Confidence level (default: 0.95)
    seed : int or None
        Random seed for reproducibility
    
    Returns
    -------
    BootstrapResult
        Point estimate and CI bounds
    """
    x = _validate_series(x, min_length=64)
    n = len(x)
    
    rng = np.random.default_rng(seed if seed is not None else 42)
    
    if block_size is None:
        block_size = max(8, int(np.sqrt(n)))
    
    n_blocks = int(np.ceil(n / block_size))
    
    # Point estimate
    point = float(estimator(x))
    
    # Bootstrap replications
    bootstrap_values = []
    for _ in range(n_bootstrap):
        # Sample block starting positions with replacement
        starts = rng.integers(0, n - block_size + 1, size=n_blocks)
        # Concatenate blocks
        boot_sample = np.concatenate([x[s:s + block_size] for s in starts])[:n]
        try:
            bootstrap_values.append(float(estimator(boot_sample)))
        except Exception:
            continue
    
    bootstrap_values = np.array(bootstrap_values)
    
    alpha = (1 - confidence) / 2
    ci_lower = float(np.quantile(bootstrap_values, alpha))
    ci_upper = float(np.quantile(bootstrap_values, 1 - alpha))
    
    return BootstrapResult(
        point_estimate=point,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        confidence_level=confidence,
        n_bootstrap=len(bootstrap_values),
        block_size=block_size,
        bootstrap_distribution=bootstrap_values
    )


# =============================================================================
# Phase-Randomized Surrogate Test
# =============================================================================

def generate_phase_surrogate(
    x: np.ndarray,
    rng: Optional[np.random.Generator] = None
) -> np.ndarray:
    """
    Generate phase-randomized surrogate preserving power spectrum.
    
    The surrogate has the same autocorrelation structure (linear
    correlations) as the original but random phases, destroying
    any nonlinear structure.
    
    Parameters
    ----------
    x : np.ndarray
        Input time series
    rng : Generator or None
        Random number generator
    
    Returns
    -------
    np.ndarray
        Phase-randomized surrogate series
    """
    x = _clean_series(x)
    n = len(x)
    
    if rng is None:
        rng = np.random.default_rng()
    
    # FFT
    X = rfft(x)
    magnitude = np.abs(X)
    phase = np.angle(X)
    
    # Randomize phases (preserving DC and Nyquist)
    random_phase = rng.uniform(-np.pi, np.pi, size=len(phase))
    random_phase[0] = phase[0]  # Keep DC component
    if len(random_phase) > 1:
        random_phase[-1] = phase[-1]  # Keep Nyquist if present
    
    # Reconstruct with random phases
    X_surrogate = magnitude * np.exp(1j * random_phase)
    surrogate = irfft(X_surrogate, n=n)
    
    return surrogate


def surrogate_test(
    x: np.ndarray,
    estimator: Callable[[np.ndarray], float],
    n_surrogates: int = 500,
    seed: Optional[int] = None
) -> SurrogateResult:
    """
    Phase-randomized surrogate test for time series estimator.
    
    Tests whether the observed estimator value differs significantly
    from what would be expected under the null hypothesis that only
    linear correlations (power spectrum) matter.
    
    Parameters
    ----------
    x : np.ndarray
        Input time series
    estimator : callable
        Function that takes array and returns scalar estimate
    n_surrogates : int
        Number of surrogate replications
    seed : int or None
        Random seed for reproducibility
    
    Returns
    -------
    SurrogateResult
        Observed value, surrogate distribution statistics, p-value
    """
    x = _validate_series(x)
    
    rng = np.random.default_rng(seed if seed is not None else 2025)
    
    # Observed value
    observed = float(estimator(x))
    
    # Surrogate distribution
    surrogate_values = []
    for _ in range(n_surrogates):
        surr = generate_phase_surrogate(x, rng)
        try:
            surrogate_values.append(float(estimator(surr)))
        except Exception:
            continue
    
    surrogate_values = np.array(surrogate_values)
    
    if len(surrogate_values) == 0:
        raise ValueError("All surrogate estimations failed")
    
    mean_surr = float(np.mean(surrogate_values))
    std_surr = float(np.std(surrogate_values)) + 1e-12  # Avoid division by zero
    
    # Two-sided p-value (rank-based)
    n_extreme = np.sum(
        np.abs(surrogate_values - mean_surr) >= np.abs(observed - mean_surr)
    )
    p_value = float((n_extreme + 1) / (len(surrogate_values) + 1))
    
    return SurrogateResult(
        observed_value=observed,
        surrogate_mean=mean_surr,
        surrogate_std=std_surr,
        p_value_two_sided=p_value,
        n_surrogates=len(surrogate_values),
        surrogate_distribution=surrogate_values
    )


# =============================================================================
# High-Level Analysis Functions
# =============================================================================

def estimate_hurst_dfa2(
    x: np.ndarray,
    fit_range: Optional[Tuple[int, int]] = None
) -> Dict[str, float]:
    """
    Convenience wrapper for DFA-2 Hurst exponent estimation.
    
    For stationary increment processes, H ≈ α (DFA scaling exponent).
    
    Returns
    -------
    dict
        Keys: alpha, H, r2
    """
    result = dfa(x, order=2, fit_range=fit_range)
    return {
        "alpha": result.alpha,
        "H": result.alpha,  # For fGn-like processes
        "r2": result.r2
    }


def full_lrd_analysis(
    x: np.ndarray,
    domain: str,
    series_name: str,
    n_bootstrap: int = 500,
    n_surrogates: int = 500,
    bootstrap_seed: int = 42,
    surrogate_seed: int = 2025
) -> LRDAnalysisResult:
    """
    Complete LRD analysis: DFA-2 + Bootstrap CI + Surrogate Test.
    
    Parameters
    ----------
    x : np.ndarray
        Input time series
    domain : str
        Domain name (e.g., "crypto", "earthquakes")
    series_name : str
        Series description
    n_bootstrap : int
        Number of bootstrap replications
    n_surrogates : int
        Number of surrogate replications
    bootstrap_seed : int
        Seed for bootstrap
    surrogate_seed : int
        Seed for surrogate test
    
    Returns
    -------
    LRDAnalysisResult
        Complete analysis results
    """
    x = _validate_series(x)
    
    # Estimator function for bootstrap and surrogate
    def hurst_estimator(series: np.ndarray) -> float:
        try:
            return dfa(series, order=2).alpha
        except Exception:
            return np.nan
    
    # DFA-2
    dfa_result = dfa(x, order=2)
    
    # Bootstrap CI
    boot_result = block_bootstrap_ci(
        x, hurst_estimator,
        n_bootstrap=n_bootstrap,
        seed=bootstrap_seed
    )
    
    # Surrogate test
    surr_result = surrogate_test(
        x, hurst_estimator,
        n_surrogates=n_surrogates,
        seed=surrogate_seed
    )
    
    # Significance flags
    is_significant = surr_result.p_value_two_sided < 0.05
    # For volatility/LRD detection: CI should not include 0.5
    is_stable = not (boot_result.ci_lower <= 0.5 <= boot_result.ci_upper)
    
    return LRDAnalysisResult(
        domain=domain,
        series_name=series_name,
        n=len(x),
        alpha=dfa_result.alpha,
        H=dfa_result.alpha,
        r2=dfa_result.r2,
        ci_lower=boot_result.ci_lower,
        ci_upper=boot_result.ci_upper,
        surrogate_mean=surr_result.surrogate_mean,
        surrogate_std=surr_result.surrogate_std,
        surrogate_p=surr_result.p_value_two_sided,
        is_significant=is_significant,
        is_stable=is_stable
    )


# =============================================================================
# Robustness Checks
# =============================================================================

def scale_range_sensitivity(
    x: np.ndarray,
    ranges: Optional[List[Tuple[int, int]]] = None
) -> Dict[str, float]:
    """
    Test sensitivity of DFA exponent to fit range selection.
    
    A stable estimate should not vary dramatically across ranges.
    """
    x = _validate_series(x)
    
    # Default: test 3 different ranges
    if ranges is None:
        n_scales = 30
        ranges = [
            (int(0.1 * n_scales), int(0.7 * n_scales)),
            (int(0.2 * n_scales), int(0.8 * n_scales)),
            (int(0.3 * n_scales), int(0.9 * n_scales)),
        ]
    
    alphas = []
    for fit_range in ranges:
        try:
            result = dfa(x, order=2, fit_range=fit_range)
            alphas.append(result.alpha)
        except Exception:
            continue
    
    if len(alphas) == 0:
        return {"mean_alpha": np.nan, "std_alpha": np.nan, "n_valid": 0}
    
    return {
        "mean_alpha": float(np.mean(alphas)),
        "std_alpha": float(np.std(alphas)),
        "range_min": float(np.min(alphas)),
        "range_max": float(np.max(alphas)),
        "n_valid": len(alphas)
    }


def rolling_dfa(
    x: np.ndarray,
    window_size: int,
    step: int = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Rolling DFA for detecting regime changes.
    
    Returns
    -------
    positions : np.ndarray
        Center positions of each window
    alphas : np.ndarray
        DFA exponents for each window
    """
    x = _validate_series(x)
    n = len(x)
    
    if step is None:
        step = window_size // 4
    
    positions = []
    alphas = []
    
    for start in range(0, n - window_size + 1, step):
        end = start + window_size
        window = x[start:end]
        try:
            result = dfa(window, order=2, min_scale=8)
            positions.append(start + window_size // 2)
            alphas.append(result.alpha)
        except Exception:
            continue
    
    return np.array(positions), np.array(alphas)
