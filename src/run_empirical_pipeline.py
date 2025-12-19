#!/usr/bin/env python3
"""
LRD v6 — Reproducible Empirical Evidence Pipeline

This script:
1. Downloads/loads real-world datasets (where permitted)
2. Preprocesses them into analysis-ready time series
3. Runs DFA-2 + bootstrap CI + surrogate tests
4. Saves results as JSON/CSV and DFA plots

Usage:
    python src/run_empirical_pipeline.py --all
    python src/run_empirical_pipeline.py --domain crypto
    python src/run_empirical_pipeline.py --domain earthquakes
    python src/run_empirical_pipeline.py --domain genomics
    python src/run_empirical_pipeline.py --domain hrv

Author: Igor Chechelnitsky
License: CC BY 4.0
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import asdict
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests
from tqdm import tqdm

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from lrd_metrics import (
    dfa, full_lrd_analysis, LRDAnalysisResult,
    scale_range_sensitivity, rolling_dfa
)


# =============================================================================
# Configuration
# =============================================================================

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"

RESULTS.mkdir(exist_ok=True, parents=True)
FIGURES.mkdir(exist_ok=True, parents=True)


# =============================================================================
# Plotting
# =============================================================================

def save_dfa_plot(
    domain: str,
    name: str,
    x: np.ndarray,
    output_path: Path,
    alpha: float = None,
    r2: float = None
) -> None:
    """Save DFA log-log scaling plot."""
    try:
        result = dfa(x, order=2)
    except Exception as e:
        print(f"  [WARN] Could not generate DFA plot: {e}")
        return
    
    s = result.scales
    F = result.fluct
    i0, i1 = result.fit_range
    
    if alpha is None:
        alpha = result.alpha
    if r2 is None:
        r2 = result.r2
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # All points
    ax.loglog(s, F, 'o-', color='steelblue', alpha=0.6, markersize=4, label='F(s)')
    
    # Fit range highlighted
    ax.loglog(s[i0:i1], F[i0:i1], 'o-', color='darkblue', markersize=6, label='Fit range')
    
    # Fit line
    fit_x = s[i0:i1]
    fit_y = np.exp(result.intercept) * fit_x ** alpha
    ax.loglog(fit_x, fit_y, '--', color='red', linewidth=2, 
              label=f'α = {alpha:.3f}, R² = {r2:.3f}')
    
    ax.set_xlabel('Scale s (window size)', fontsize=12)
    ax.set_ylabel('F(s) — Fluctuation', fontsize=12)
    ax.set_title(f'{domain}: DFA-2 Scaling ({name})', fontsize=14)
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_path}")


def save_rolling_dfa_plot(
    domain: str,
    name: str,
    positions: np.ndarray,
    alphas: np.ndarray,
    output_path: Path
) -> None:
    """Save rolling DFA exponent plot for regime analysis."""
    fig, ax = plt.subplots(figsize=(10, 4))
    
    ax.plot(positions, alphas, 'o-', color='steelblue', markersize=3)
    ax.axhline(y=0.5, color='gray', linestyle='--', label='H = 0.5 (memoryless)')
    ax.axhline(y=np.mean(alphas), color='red', linestyle='-', alpha=0.7,
               label=f'Mean α = {np.mean(alphas):.3f}')
    
    ax.set_xlabel('Position (data point index)', fontsize=12)
    ax.set_ylabel('DFA Exponent α', fontsize=12)
    ax.set_title(f'{domain}: Rolling DFA-2 ({name})', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_path}")


# =============================================================================
# Domain Loaders
# =============================================================================

def load_crypto_yahoo(
    symbol: str = "BTC-USD",
    start: str = "2016-01-01",
    end: str = "2024-12-01"
) -> Dict[str, np.ndarray]:
    """
    Download OHLCV from Yahoo Finance.
    
    Returns:
        returns: log returns
        volproxy: |log returns| (volatility proxy)
        prices: adjusted close prices
    """
    print(f"  Downloading {symbol} from Yahoo Finance...")
    
    import time
    from datetime import datetime as dt
    
    def to_unix(s: str) -> int:
        return int(time.mktime(dt.strptime(s, "%Y-%m-%d").timetuple()))
    
    p1 = to_unix(start)
    p2 = to_unix(end)
    
    url = f"https://query1.finance.yahoo.com/v7/finance/download/{symbol}"
    params = {
        "period1": p1,
        "period2": p2,
        "interval": "1d",
        "events": "history",
        "includeAdjustedClose": "true"
    }
    
    headers = {
        "User-Agent": "Mozilla/5.0 (compatible; LRD-v6/1.0)"
    }
    
    r = requests.get(url, params=params, headers=headers, timeout=60)
    r.raise_for_status()
    
    df = pd.read_csv(StringIO(r.text))
    df = df.dropna(subset=["Adj Close"])
    
    price = df["Adj Close"].values.astype(float)
    
    if len(price) < 100:
        raise ValueError(f"Insufficient data: only {len(price)} points")
    
    # Log returns
    log_returns = np.diff(np.log(price))
    
    # Volatility proxy: absolute log returns
    vol_proxy = np.abs(log_returns)
    
    print(f"  Downloaded {len(price)} data points ({start} to {end})")
    
    return {
        "returns": log_returns,
        "volproxy": vol_proxy,
        "prices": price
    }


def load_earthquakes_usgs(
    start: str = "2010-01-01",
    end: str = "2023-01-01",
    minmag: float = 4.5,
    limit: int = 20000
) -> np.ndarray:
    """
    Download earthquake catalog from USGS and return daily event counts.
    
    Source: USGS Earthquake Hazards Program (FDSN Event Web Service)
    """
    print(f"  Downloading USGS earthquake catalog (M≥{minmag})...")
    
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query"
    params = {
        "format": "geojson",
        "starttime": start,
        "endtime": end,
        "minmagnitude": minmag,
        "orderby": "time-asc",
        "limit": limit
    }
    
    r = requests.get(url, params=params, timeout=120)
    r.raise_for_status()
    data = r.json()
    
    times = []
    for feat in data.get("features", []):
        ts = feat["properties"].get("time")
        if ts is None:
            continue
        times.append(pd.to_datetime(int(ts), unit="ms", utc=True).date())
    
    if len(times) < 200:
        raise ValueError(f"Too few events: {len(times)}. Adjust parameters.")
    
    # Aggregate to daily counts
    series = pd.Series(1, index=pd.to_datetime(times))
    daily = series.resample("D").sum().fillna(0.0)
    x = daily.values.astype(float)
    
    print(f"  Downloaded {len(times)} events → {len(x)} daily counts")
    
    return x


def load_genomics_ncbi(
    accession: str = "NC_000001.11",
    max_len: int = 200000
) -> np.ndarray:
    """
    Download FASTA sequence from NCBI and encode as +1/-1 series.
    
    Encoding:
        Purine (A, G) → +1
        Pyrimidine (C, T) → -1
    """
    print(f"  Downloading NCBI sequence {accession}...")
    
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        raise ImportError("Biopython required: pip install biopython")
    
    Entrez.email = "lrd-v6@example.com"  # Required by NCBI
    
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text"
    )
    record = SeqIO.read(handle, "fasta")
    handle.close()
    
    seq = str(record.seq).upper()
    
    if max_len is not None:
        seq = seq[:max_len]
    
    # Encode: purine +1, pyrimidine -1
    mapping = {"A": 1, "G": 1, "C": -1, "T": -1}
    x = np.array([mapping.get(ch, 0) for ch in seq], dtype=float)
    x = x[x != 0]  # Remove unknowns (N, etc.)
    
    # Center the series
    x = x - np.mean(x)
    
    print(f"  Downloaded {len(seq)} bases → {len(x)} encoded points")
    
    return x


def load_hrv_physionet_example() -> np.ndarray:
    """
    Load example HRV (RR intervals) from PhysioNet.
    
    NOTE: This is a stub that generates synthetic data for demonstration.
    For real analysis, use wfdb to load actual PhysioNet records.
    """
    print("  Generating synthetic HRV example (for demonstration)...")
    print("  NOTE: For real analysis, use wfdb to load PhysioNet records")
    
    # Generate synthetic RR intervals with LRD-like structure
    # This is for demonstration only - use real data for publication
    np.random.seed(42)
    n = 4096
    
    # Simple fGn-like approximation via spectral method
    freqs = np.fft.rfftfreq(n)
    freqs[0] = 1e-10  # Avoid division by zero
    
    H = 0.75  # Typical healthy HRV exponent
    power = freqs ** (-2 * H - 1)
    
    phases = np.random.uniform(-np.pi, np.pi, len(power))
    X = np.sqrt(power) * np.exp(1j * phases)
    X[0] = 0  # Zero mean
    
    rr = np.fft.irfft(X, n=n)
    rr = rr - np.mean(rr)
    
    print(f"  Generated {len(rr)} synthetic RR intervals (H ≈ {H})")
    
    return rr


# =============================================================================
# Domain Runners
# =============================================================================

def run_crypto(verbose: bool = True) -> List[LRDAnalysisResult]:
    """Analyze Bitcoin/crypto data."""
    results = []
    
    if verbose:
        print("\n" + "="*60)
        print("DOMAIN: Crypto/Finance (Bitcoin)")
        print("="*60)
    
    try:
        data = load_crypto_yahoo(symbol="BTC-USD", start="2016-01-01", end="2024-12-01")
    except Exception as e:
        print(f"  [ERROR] Failed to load crypto data: {e}")
        return results
    
    # 1. Log returns (expected: H ≈ 0.5)
    if verbose:
        print("\n  Analyzing: BTC log-returns...")
    
    try:
        result_returns = full_lrd_analysis(
            data["returns"],
            domain="crypto",
            series_name="BTC log-returns",
            n_bootstrap=500,
            n_surrogates=500
        )
        results.append(result_returns)
        
        if verbose:
            print(f"    α = {result_returns.alpha:.4f}, H = {result_returns.H:.4f}")
            print(f"    CI: [{result_returns.ci_lower:.4f}, {result_returns.ci_upper:.4f}]")
            print(f"    Surrogate p = {result_returns.surrogate_p:.4f}")
        
        save_dfa_plot(
            "crypto", "BTC log-returns",
            data["returns"],
            FIGURES / "crypto_btc_returns_dfa.png"
        )
    except Exception as e:
        print(f"    [ERROR] Returns analysis failed: {e}")
    
    # 2. Volatility proxy (expected: H > 0.5)
    if verbose:
        print("\n  Analyzing: BTC |log-returns| (volatility proxy)...")
    
    try:
        result_vol = full_lrd_analysis(
            data["volproxy"],
            domain="crypto",
            series_name="BTC |log-returns| (volatility proxy)",
            n_bootstrap=500,
            n_surrogates=500
        )
        results.append(result_vol)
        
        if verbose:
            print(f"    α = {result_vol.alpha:.4f}, H = {result_vol.H:.4f}")
            print(f"    CI: [{result_vol.ci_lower:.4f}, {result_vol.ci_upper:.4f}]")
            print(f"    Surrogate p = {result_vol.surrogate_p:.4f}")
        
        save_dfa_plot(
            "crypto", "BTC volatility proxy",
            data["volproxy"],
            FIGURES / "crypto_btc_volproxy_dfa.png"
        )
    except Exception as e:
        print(f"    [ERROR] Volatility analysis failed: {e}")
    
    # 3. Rolling DFA for regime detection (on volatility)
    if verbose:
        print("\n  Computing rolling DFA (volatility, 500-day windows)...")
    
    try:
        if len(data["volproxy"]) > 1000:
            positions, alphas = rolling_dfa(data["volproxy"], window_size=500, step=50)
            save_rolling_dfa_plot(
                "crypto", "BTC volatility",
                positions, alphas,
                FIGURES / "crypto_btc_rolling_dfa.png"
            )
    except Exception as e:
        print(f"    [ERROR] Rolling DFA failed: {e}")
    
    return results


def run_earthquakes(verbose: bool = True) -> List[LRDAnalysisResult]:
    """Analyze earthquake data."""
    results = []
    
    if verbose:
        print("\n" + "="*60)
        print("DOMAIN: Geophysics (Earthquakes)")
        print("="*60)
    
    try:
        x = load_earthquakes_usgs(start="2010-01-01", end="2023-01-01", minmag=4.5)
    except Exception as e:
        print(f"  [ERROR] Failed to load earthquake data: {e}")
        return results
    
    if verbose:
        print("\n  Analyzing: USGS daily earthquake counts (M≥4.5)...")
    
    try:
        result = full_lrd_analysis(
            x,
            domain="earthquakes",
            series_name="USGS daily counts (M≥4.5)",
            n_bootstrap=500,
            n_surrogates=500
        )
        results.append(result)
        
        if verbose:
            print(f"    α = {result.alpha:.4f}, H = {result.H:.4f}")
            print(f"    CI: [{result.ci_lower:.4f}, {result.ci_upper:.4f}]")
            print(f"    Surrogate p = {result.surrogate_p:.4f}")
        
        save_dfa_plot(
            "earthquakes", "USGS daily counts",
            x,
            FIGURES / "earthquakes_usgs_dfa.png"
        )
    except Exception as e:
        print(f"    [ERROR] Analysis failed: {e}")
    
    return results


def run_genomics(verbose: bool = True) -> List[LRDAnalysisResult]:
    """Analyze genomics data."""
    results = []
    
    if verbose:
        print("\n" + "="*60)
        print("DOMAIN: Genomics (DNA sequences)")
        print("="*60)
    
    try:
        x = load_genomics_ncbi(accession="NC_000001.11", max_len=200000)
    except Exception as e:
        print(f"  [ERROR] Failed to load genomics data: {e}")
        return results
    
    if verbose:
        print("\n  Analyzing: Purine/pyrimidine walk (+1/-1)...")
    
    try:
        result = full_lrd_analysis(
            x,
            domain="genomics",
            series_name="NCBI NC_000001.11 purine/pyrimidine (+1/-1)",
            n_bootstrap=500,
            n_surrogates=500
        )
        results.append(result)
        
        if verbose:
            print(f"    α = {result.alpha:.4f}, H = {result.H:.4f}")
            print(f"    CI: [{result.ci_lower:.4f}, {result.ci_upper:.4f}]")
            print(f"    Surrogate p = {result.surrogate_p:.4f}")
        
        save_dfa_plot(
            "genomics", "DNA purine/pyrimidine walk",
            x,
            FIGURES / "genomics_dna_dfa.png"
        )
    except Exception as e:
        print(f"    [ERROR] Analysis failed: {e}")
    
    return results


def run_hrv(verbose: bool = True) -> List[LRDAnalysisResult]:
    """Analyze HRV data."""
    results = []
    
    if verbose:
        print("\n" + "="*60)
        print("DOMAIN: Physiology (Heart Rate Variability)")
        print("="*60)
    
    try:
        x = load_hrv_physionet_example()
    except Exception as e:
        print(f"  [ERROR] Failed to load HRV data: {e}")
        return results
    
    if verbose:
        print("\n  Analyzing: RR interval deviations...")
    
    try:
        result = full_lrd_analysis(
            x,
            domain="hrv",
            series_name="RR intervals (synthetic demo)",
            n_bootstrap=500,
            n_surrogates=500
        )
        results.append(result)
        
        if verbose:
            print(f"    α = {result.alpha:.4f}, H = {result.H:.4f}")
            print(f"    CI: [{result.ci_lower:.4f}, {result.ci_upper:.4f}]")
            print(f"    Surrogate p = {result.surrogate_p:.4f}")
        
        save_dfa_plot(
            "hrv", "RR intervals",
            x,
            FIGURES / "hrv_rr_dfa.png"
        )
    except Exception as e:
        print(f"    [ERROR] Analysis failed: {e}")
    
    return results


# =============================================================================
# Output Writers
# =============================================================================

def write_results(results: List[LRDAnalysisResult]) -> Tuple[Path, Path]:
    """Write results to JSON and CSV."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # JSON
    json_path = RESULTS / "metrics_summary.json"
    output = {
        "generated": timestamp,
        "package": "LRD v6 — Empirical Evidence",
        "author": "Igor Chechelnitsky",
        "n_results": len(results),
        "results": [r.to_dict() for r in results]
    }
    json_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    
    # CSV
    csv_path = RESULTS / "metrics_table.csv"
    df = pd.DataFrame([r.to_dict() for r in results])
    df.to_csv(csv_path, index=False)
    
    return json_path, csv_path


def write_provenance(domains_run: List[str]) -> Path:
    """Write provenance metadata."""
    prov_path = RESULTS / "provenance.json"
    
    provenance = {
        "generated": datetime.now().isoformat(),
        "package_version": "6.0.0",
        "domains_analyzed": domains_run,
        "data_sources": {
            "crypto": {
                "source": "Yahoo Finance",
                "symbol": "BTC-USD",
                "endpoint": "query1.finance.yahoo.com",
                "license": "Yahoo Terms of Service"
            },
            "earthquakes": {
                "source": "USGS Earthquake Hazards Program",
                "endpoint": "earthquake.usgs.gov/fdsnws/event/1/query",
                "format": "GeoJSON",
                "license": "Public Domain"
            },
            "genomics": {
                "source": "NCBI GenBank",
                "accession": "NC_000001.11",
                "endpoint": "NCBI Entrez",
                "license": "Public Domain"
            },
            "hrv": {
                "source": "Synthetic (demonstration)",
                "note": "Use wfdb + PhysioNet for real data",
                "recommended": "PhysioNet NSR2DB or CHF2DB"
            }
        },
        "methodology": {
            "estimator": "DFA-2 (quadratic detrending)",
            "uncertainty": "Block bootstrap (500 replications)",
            "hypothesis_test": "Phase-randomized surrogates (500 replications)",
            "reference": "Peng et al. (1994); Kantelhardt et al. (2002)"
        }
    }
    
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    
    return prov_path


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="LRD v6 — Empirical Evidence Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python src/run_empirical_pipeline.py --all
  python src/run_empirical_pipeline.py --domain crypto
  python src/run_empirical_pipeline.py --domain earthquakes --domain genomics
        """
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Run all domains"
    )
    parser.add_argument(
        "--domain", action="append", dest="domains",
        choices=["crypto", "earthquakes", "genomics", "hrv"],
        help="Domain to analyze (can be specified multiple times)"
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress output"
    )
    
    args = parser.parse_args()
    
    if args.all:
        domains = ["crypto", "earthquakes", "genomics", "hrv"]
    elif args.domains:
        domains = args.domains
    else:
        parser.error("Provide --all or at least one --domain")
    
    verbose = not args.quiet
    
    if verbose:
        print("="*60)
        print("LRD v6 — Empirical Evidence Pipeline")
        print("="*60)
        print(f"Domains: {', '.join(domains)}")
        print(f"Output: {RESULTS}")
    
    # Run domains
    all_results: List[LRDAnalysisResult] = []
    domains_run = []
    
    domain_runners = {
        "crypto": run_crypto,
        "earthquakes": run_earthquakes,
        "genomics": run_genomics,
        "hrv": run_hrv
    }
    
    for domain in domains:
        runner = domain_runners[domain]
        try:
            results = runner(verbose=verbose)
            all_results.extend(results)
            domains_run.append(domain)
        except Exception as e:
            print(f"\n[ERROR] Domain '{domain}' failed: {e}")
    
    if not all_results:
        print("\n[ERROR] No results produced. Check errors above.")
        sys.exit(1)
    
    # Write outputs
    if verbose:
        print("\n" + "="*60)
        print("SAVING RESULTS")
        print("="*60)
    
    json_path, csv_path = write_results(all_results)
    prov_path = write_provenance(domains_run)
    
    print(f"\nSaved: {json_path}")
    print(f"Saved: {csv_path}")
    print(f"Saved: {prov_path}")
    print(f"Figures: {FIGURES}")
    
    # Summary
    if verbose:
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        print(f"{'Domain':<15} {'Series':<40} {'H':>8} {'CI':>20} {'p':>8}")
        print("-"*95)
        for r in all_results:
            ci_str = f"[{r.ci_lower:.3f}, {r.ci_upper:.3f}]"
            print(f"{r.domain:<15} {r.series_name[:38]:<40} {r.H:>8.4f} {ci_str:>20} {r.surrogate_p:>8.4f}")
    
    print("\nDone.")


if __name__ == "__main__":
    main()
