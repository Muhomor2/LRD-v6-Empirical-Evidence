# LRD v6 — Empirical Evidence for Long-Range Dependence / Fractal Memory

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Companion to**: [LRD v5 – Monte Carlo Validation](https://zenodo.org/records/17800770) (DOI: 10.5281/zenodo.17800770)

**Author**: Igor Chechelnitsky — Independent Researcher, Israel (Ashkelon)  
**ORCID**: [0009-0007-4607-1946](https://orcid.org/0009-0007-4607-1946)

---

## Abstract

This package extends **LRD v5** (Monte-Carlo validation on synthetic fGn) with a **reproducible, multi-domain empirical validation** of long-range dependence (LRD) / fractal memory using **real-world datasets**.

**Central Claim (Falsifiable)**: There exist real-world time series for which a conservative DFA-2 pipeline yields scaling exponents consistent with LRD beyond what is expected from (a) simple trends, (b) short-memory ARMA structure, or (c) preserved-spectrum surrogates.

**Expected Pattern (Finance/Bitcoin)**:
- **Returns**: Often near memoryless (H ≈ 0.5)
- **Volatility proxy |r_t|**: Persistent dependence (H > 0.5), consistent with volatility clustering

This asymmetry between returns and volatility is the canonical empirical signature in financial markets.

---

## Domains Covered

| Domain | Dataset | Series Type | Expected H |
|--------|---------|-------------|------------|
| **Finance/Crypto** | BTC-USD (Yahoo Finance) | Log-returns | ~0.5 (weak memory) |
| **Finance/Crypto** | BTC-USD (Yahoo Finance) | \|log-returns\| | >0.5 (volatility clustering) |
| **Geophysics** | USGS Earthquake Catalog | Daily event counts (M≥4.5) | >0.5 (aftershock clustering) |
| **Physiology** | PhysioNet RR Intervals | RR deviation series | Variable (healthy vs pathology) |
| **Genomics** | NCBI GenBank FASTA | Purine/pyrimidine (+1/-1) | ~0.6–0.8 (sequence correlations) |

---

## Methodology

### DFA-2 (Detrended Fluctuation Analysis, quadratic detrending)

1. Convert to analysis series (returns, RR deviations, daily counts, ±1 mapping)
2. Build integrated profile: Y(k) = Σᵢ₌₁ᵏ (xᵢ - ⟨x⟩)
3. For scale s, detrend segments with polynomial of order 2, compute RMS fluctuation F(s)
4. Fit log F(s) vs log s on pre-registered scale range; slope = α
5. Interpret α as scaling exponent; map to H under appropriate assumptions

### Uncertainty Quantification

- **Block Bootstrap CI**: Preserves temporal dependence within blocks (95% CI)
- **Phase-Randomized Surrogate Test**: Destroys nonlinear/phase structure while preserving power spectrum

### Robustness Controls

- Scale-range sensitivity (3+ ranges)
- Rolling windows for Bitcoin (regime variation)
- Returns vs volatility separation
- Time-reversed controls

---

## Falsification Criteria

Evidence is **not accepted** if:

1. Surrogate tests preserve the scaling exponent (p > 0.05)
2. Bootstrap CI includes null region (H = 0.5 for memoryless)
3. Time-reversed controls show different behavior
4. Scale-range perturbations cause instability

---

## Quick Start

### 1) Create Python Environment

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### 2) Run All Domains

```bash
python src/run_empirical_pipeline.py --all
```

### 3) Run Single Domain

```bash
python src/run_empirical_pipeline.py --domain crypto
python src/run_empirical_pipeline.py --domain earthquakes
python src/run_empirical_pipeline.py --domain genomics
python src/run_empirical_pipeline.py --domain hrv
```

### Outputs

- `results/metrics_summary.json` — Machine-readable results
- `results/metrics_table.csv` — Tabular summary
- `results/figures/*.png` — DFA log-log scaling plots

---

## AI Relevance

### Forecasting in Markets with Memory

If volatility proxies exhibit persistent dependence, risk is not memoryless. Any AI system forecasting risk, stress, or tail exposure must treat volatility as a stateful process rather than white noise.

### Agentic AI and Feedback Loops

Automated trading systems can induce reflexive dynamics: strategies adapt to volatility states; volatility states influence strategy selection; the system can amplify clustering.

### Robustness Benchmark

**Proposed benchmark**: Any AI market model claiming realism should reproduce:
- Returns: weak/unstable memory (H ≈ 0.5)
- Volatility: strong/stable memory (H > 0.5)

Failure to reproduce this is a red flag for realism.

---

## Data Licensing

This package **does not redistribute raw third-party datasets** by default. We store only:
- Derived metrics (H/α, CI, surrogate p-values)
- Plots (DFA log-log scaling)
- Provenance metadata (query parameters, accessions, timestamps)

See `docs/DATA_SOURCES.md` for dataset sources and `docs/THIRD_PARTY_NOTICES.md` for licensing notes.

---

## Citation

If you use this software or derived results, please cite:

```bibtex
@software{chechelnitsky_lrd_v6_2025,
  author       = {Chechelnitsky, Igor},
  title        = {{LRD v6: Empirical Evidence for Long-Range Dependence / Fractal Memory}},
  year         = 2025,
  publisher    = {Zenodo},
  version      = {6.0.0},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

Also cite: [LRD v5 – Monte Carlo Validation](https://zenodo.org/records/17800770) (DOI: 10.5281/zenodo.17800770)

---

## Repository Structure

```
LRD-v6-Empirical-Evidence/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── LICENSE                      # CC BY 4.0
├── CITATION.cff                 # Citation metadata
├── zenodo.json                  # Zenodo metadata
├── src/
│   ├── lrd_metrics.py          # Core DFA-2, bootstrap, surrogate tests
│   └── run_empirical_pipeline.py  # Main pipeline script
├── docs/
│   ├── DATA_SOURCES.md         # Dataset sources and provenance
│   ├── THIRD_PARTY_NOTICES.md  # Third-party licensing
│   ├── AUDIT_CHECKLIST.md      # Reviewer audit checklist
│   └── LICENSE-OSL-ER-v1.txt   # Ethical research policy
├── latex/
│   └── LRD_v6_Empirical_Evidence.tex  # Academic manuscript
├── tests/
│   └── test_lrd_metrics.py     # Unit tests
└── results/                     # Generated outputs (gitignored)
    ├── metrics_summary.json
    ├── metrics_table.csv
    └── figures/
```

---

## References

1. Peng C-K et al. (1994). Quantification of scaling exponents and crossover phenomena in nonstationary heartbeat time series. Chaos.
2. Ding Z, Granger CWJ, Engle RF (1993). A long memory property of stock market returns and a new model. J Empirical Finance.
3. Bariviera AF (2017). The inefficiency of Bitcoin revisited: A dynamic approach. Economics Letters.
4. Kantelhardt JW et al. (2002). Multifractal detrended fluctuation analysis of nonstationary time series. Physica A.

---

## License

This work is licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

Third-party datasets retain their original licenses — see `docs/THIRD_PARTY_NOTICES.md`.
