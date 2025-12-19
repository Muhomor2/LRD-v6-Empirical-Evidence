# Audit Checklist — LRD v6 Empirical Evidence

This checklist is for reviewers to verify the scientific validity and reproducibility of the empirical evidence presented.

---

## ✅ Minimal Acceptance Criteria

For results to be accepted as valid evidence of long-range dependence:

### Methodology

- [ ] **DFA-2 exponent (α) reported** with fit quality (R²)
- [ ] **Block-bootstrap CI reported** (95% confidence interval)
- [ ] **Phase-randomized surrogate test** p-value reported
- [ ] **Scale-range sensitivity** tested (at least 3 ranges)

### Finance/Crypto Specific

- [ ] **Returns vs volatility reported separately**
  - Returns: expected H ≈ 0.5 (weak/no memory)
  - Volatility: expected H > 0.5 (clustering)
- [ ] **Rolling windows shown** (to avoid cherry-picking)

### Reproducibility

- [ ] **Dataset provenance documented** (sources, query parameters)
- [ ] **License notes included** for all data sources
- [ ] **Code is runnable** with documented dependencies
- [ ] **Random seeds specified** for bootstrap/surrogate

---

## ⚠️ Key Failure Modes to Check

### 1. Nonstationarity Masquerading as LRD

**Risk**: Trends in raw data inflate H estimates

**Checks**:
- [ ] Series appropriately transformed (returns, differences, centering)
- [ ] DFA-2 detrending applied (quadratic, not just linear)
- [ ] Rolling windows show consistent estimates

### 2. Scale-Range Cherry Picking

**Risk**: Fitting only the "good" part of the log-log plot

**Checks**:
- [ ] Multiple scale ranges tested
- [ ] Results stable across range perturbations
- [ ] Central 60% default range is reasonable

### 3. Heavy Tails and Outliers

**Risk**: Extreme values distort fluctuation estimates

**Checks**:
- [ ] Data inspected for extreme outliers
- [ ] Sensitivity to winsorization documented (if applicable)

### 4. Microstructure Noise (High-Frequency Data)

**Risk**: Market microstructure creates spurious correlations

**Checks**:
- [ ] Daily data used for baseline (avoids intraday noise)
- [ ] If using high-frequency, noise-robust estimators applied

### 5. Regime Shifts

**Risk**: Single exponent doesn't represent time-varying behavior

**Checks**:
- [ ] Rolling DFA shows regime variation
- [ ] Not claiming "permanent" LRD structure
- [ ] Time range and regime context documented

### 6. Surrogate Test Survival

**Risk**: LRD is trivially spectral (linear correlations only)

**Checks**:
- [ ] Observed exponent differs from surrogate distribution (p < 0.05)
- [ ] If p > 0.05, evidence is weak and noted as such

---

## 📊 Expected Results Pattern

A skeptical reviewer should expect:

### Finance/Crypto (Bitcoin)

| Series | Expected H | Interpretation |
|--------|------------|----------------|
| Log returns | ~0.5 | Weak/no memory (near random walk) |
| \|log returns\| | >0.5 | Volatility clustering (long memory) |

This asymmetry is the canonical empirical signature.

### Geophysics (Earthquakes)

| Series | Expected H | Interpretation |
|--------|------------|----------------|
| Daily counts | >0.5 | Aftershock clustering |

### Physiology (HRV)

| Series | Expected H | Interpretation |
|--------|------------|----------------|
| Healthy RR | ~0.8-1.0 | Scale-invariant dynamics |
| Pathology RR | Different | Often reduced complexity |

### Genomics (DNA)

| Series | Expected H | Interpretation |
|--------|------------|----------------|
| Purine/pyrimidine | ~0.6-0.8 | Sequence correlations |

---

## 🔬 Falsification Criteria

Evidence is **NOT accepted** if:

1. ❌ Surrogate test p-value > 0.05
2. ❌ Bootstrap CI includes null region (H = 0.5 when testing for LRD)
3. ❌ Scale-range perturbations cause >0.1 change in α
4. ❌ Time-reversed series shows different behavior (for reversible processes)
5. ❌ Rolling windows show high variance (>0.15 std in α)

---

## 📝 Reviewer Notes Template

```
Domain: _______________
Series: _______________
Date Reviewed: _______________

DFA-2 Exponent (α): _______________
R²: _______________
Bootstrap CI: [___________, ___________]
Surrogate p-value: _______________

Checks Passed:
[ ] Returns vs volatility separated (if finance)
[ ] Scale-range sensitivity OK
[ ] Rolling windows consistent (if applicable)
[ ] Surrogate test significant
[ ] Provenance documented

Notes:
_______________________________________
_______________________________________

Verdict: [ ] Accept  [ ] Reject  [ ] Revise
```

---

## References

1. Peng CK et al. (1994). Quantification of scaling exponents... *Chaos*.
2. Kantelhardt JW et al. (2002). Multifractal detrended fluctuation analysis... *Physica A*.
3. Ding Z, Granger CWJ, Engle RF (1993). A long memory property of stock market returns... *J Empirical Finance*.
4. Bariviera AF (2017). The inefficiency of Bitcoin revisited... *Economics Letters*.
