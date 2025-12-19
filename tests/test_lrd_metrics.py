"""
Unit tests for LRD v6 metrics module.

Run with: pytest tests/test_lrd_metrics.py -v
"""

import numpy as np
import pytest
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lrd_metrics import (
    dfa, DFAResult,
    block_bootstrap_ci, BootstrapResult,
    surrogate_test, SurrogateResult,
    generate_phase_surrogate,
    estimate_hurst_dfa2,
    full_lrd_analysis,
    scale_range_sensitivity,
    rolling_dfa
)


class TestDFA:
    """Tests for DFA implementation."""
    
    def test_dfa_white_noise(self):
        """White noise should give H ≈ 0.5."""
        np.random.seed(42)
        x = np.random.randn(4096)
        result = dfa(x, order=2)
        
        assert isinstance(result, DFAResult)
        assert 0.4 < result.alpha < 0.6, f"Expected H ≈ 0.5, got {result.alpha}"
        assert result.r2 > 0.95, f"Poor fit quality: R² = {result.r2}"
    
    def test_dfa_random_walk(self):
        """Random walk (cumsum of white noise) should give H ≈ 1.5."""
        np.random.seed(42)
        x = np.cumsum(np.random.randn(4096))
        result = dfa(x, order=2)
        
        assert 1.3 < result.alpha < 1.7, f"Expected H ≈ 1.5, got {result.alpha}"
    
    def test_dfa_synthetic_fgn(self):
        """Synthetic fGn-like series should recover approximate H."""
        np.random.seed(42)
        n = 4096
        H_target = 0.75
        
        # Simple spectral method for fGn approximation
        freqs = np.fft.rfftfreq(n)
        freqs[0] = 1e-10
        power = freqs ** (-2 * H_target - 1)
        phases = np.random.uniform(-np.pi, np.pi, len(power))
        X = np.sqrt(power) * np.exp(1j * phases)
        X[0] = 0
        x = np.fft.irfft(X, n=n)
        
        result = dfa(x, order=2)
        
        # Allow ±0.15 tolerance for spectral approximation
        assert abs(result.alpha - H_target) < 0.15, \
            f"Expected H ≈ {H_target}, got {result.alpha}"
    
    def test_dfa_short_series_raises(self):
        """Series too short should raise ValueError."""
        x = np.random.randn(50)
        with pytest.raises(ValueError):
            dfa(x, order=2, min_scale=16)
    
    def test_dfa_all_nan_raises(self):
        """All-NaN series should raise ValueError."""
        x = np.array([np.nan] * 1000)
        with pytest.raises(ValueError):
            dfa(x, order=2)
    
    def test_dfa_handles_nan(self):
        """NaN values should be cleaned."""
        np.random.seed(42)
        x = np.random.randn(2000)
        x[100:110] = np.nan  # Insert some NaNs
        
        result = dfa(x, order=2)
        assert np.isfinite(result.alpha)


class TestBootstrap:
    """Tests for block bootstrap CI."""
    
    def test_bootstrap_basic(self):
        """Basic bootstrap should return valid CI."""
        np.random.seed(42)
        x = np.random.randn(1000)
        
        def mean_estimator(series):
            return np.mean(series)
        
        result = block_bootstrap_ci(x, mean_estimator, n_bootstrap=100, seed=42)
        
        assert isinstance(result, BootstrapResult)
        assert result.ci_lower < result.point_estimate < result.ci_upper
        assert result.confidence_level == 0.95
    
    def test_bootstrap_contains_true_value(self):
        """CI should contain true parameter with high probability."""
        np.random.seed(42)
        true_mean = 5.0
        x = np.random.randn(1000) + true_mean
        
        result = block_bootstrap_ci(x, np.mean, n_bootstrap=200, seed=42)
        
        assert result.ci_lower < true_mean < result.ci_upper
    
    def test_bootstrap_short_series_raises(self):
        """Short series should raise ValueError."""
        x = np.random.randn(30)
        with pytest.raises(ValueError):
            block_bootstrap_ci(x, np.mean)


class TestSurrogate:
    """Tests for phase-randomized surrogate tests."""
    
    def test_surrogate_preserves_spectrum(self):
        """Surrogate should have similar power spectrum."""
        np.random.seed(42)
        x = np.random.randn(1024)
        
        surr = generate_phase_surrogate(x)
        
        # Compare power spectra
        orig_power = np.abs(np.fft.rfft(x)) ** 2
        surr_power = np.abs(np.fft.rfft(surr)) ** 2
        
        # Should be approximately equal
        correlation = np.corrcoef(orig_power, surr_power)[0, 1]
        assert correlation > 0.99, f"Power spectra differ: corr = {correlation}"
    
    def test_surrogate_destroys_phase(self):
        """Surrogate should have different phase structure."""
        np.random.seed(42)
        x = np.random.randn(1024)
        
        surr = generate_phase_surrogate(x)
        
        # Series should be different
        correlation = np.corrcoef(x, surr)[0, 1]
        assert abs(correlation) < 0.3, f"Surrogate too similar: corr = {correlation}"
    
    def test_surrogate_test_basic(self):
        """Basic surrogate test should return valid result."""
        np.random.seed(42)
        x = np.random.randn(1024)
        
        result = surrogate_test(x, np.mean, n_surrogates=50, seed=42)
        
        assert isinstance(result, SurrogateResult)
        assert 0 <= result.p_value_two_sided <= 1
        assert len(result.surrogate_distribution) == 50
    
    def test_surrogate_test_white_noise(self):
        """White noise DFA exponent should not differ from surrogates."""
        np.random.seed(42)
        x = np.random.randn(2048)
        
        def hurst_est(s):
            try:
                return dfa(s, order=2).alpha
            except:
                return np.nan
        
        result = surrogate_test(x, hurst_est, n_surrogates=50, seed=42)
        
        # White noise should not be significantly different
        # (p-value should not be very small)
        assert result.p_value_two_sided > 0.01


class TestHighLevel:
    """Tests for high-level analysis functions."""
    
    def test_estimate_hurst_dfa2(self):
        """Wrapper function should return dict with correct keys."""
        np.random.seed(42)
        x = np.random.randn(2048)
        
        result = estimate_hurst_dfa2(x)
        
        assert "alpha" in result
        assert "H" in result
        assert "r2" in result
        assert result["alpha"] == result["H"]  # For fGn-like series
    
    def test_full_lrd_analysis(self):
        """Full analysis should return complete result."""
        np.random.seed(42)
        x = np.random.randn(2048)
        
        result = full_lrd_analysis(
            x, domain="test", series_name="random",
            n_bootstrap=50, n_surrogates=50
        )
        
        assert result.domain == "test"
        assert result.series_name == "random"
        assert result.n == 2048
        assert np.isfinite(result.alpha)
        assert np.isfinite(result.ci_lower)
        assert np.isfinite(result.ci_upper)
        assert np.isfinite(result.surrogate_p)
    
    def test_scale_range_sensitivity(self):
        """Scale range sensitivity should return stability metrics."""
        np.random.seed(42)
        x = np.random.randn(4096)
        
        result = scale_range_sensitivity(x)
        
        assert "mean_alpha" in result
        assert "std_alpha" in result
        assert result["n_valid"] >= 2
        
        # Stable series should have low std
        assert result["std_alpha"] < 0.1
    
    def test_rolling_dfa(self):
        """Rolling DFA should return positions and alphas."""
        np.random.seed(42)
        x = np.random.randn(4096)
        
        positions, alphas = rolling_dfa(x, window_size=512, step=128)
        
        assert len(positions) == len(alphas)
        assert len(positions) > 5
        assert all(np.isfinite(alphas))


class TestEdgeCases:
    """Tests for edge cases and robustness."""
    
    def test_constant_series_raises(self):
        """Constant series should raise or return extreme values."""
        x = np.ones(1000)
        with pytest.raises((ValueError, ZeroDivisionError)):
            dfa(x, order=2)
    
    def test_mixed_nan_inf(self):
        """Mixed NaN and Inf should be cleaned."""
        np.random.seed(42)
        x = np.random.randn(2000)
        x[100] = np.nan
        x[200] = np.inf
        x[300] = -np.inf
        
        result = dfa(x, order=2)
        assert np.isfinite(result.alpha)
    
    def test_reproducibility(self):
        """Same seed should give same results."""
        np.random.seed(42)
        x = np.random.randn(2048)
        
        result1 = block_bootstrap_ci(x, np.mean, n_bootstrap=100, seed=123)
        result2 = block_bootstrap_ci(x, np.mean, n_bootstrap=100, seed=123)
        
        assert result1.ci_lower == result2.ci_lower
        assert result1.ci_upper == result2.ci_upper


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
