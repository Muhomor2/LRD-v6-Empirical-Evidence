# Addendum (v6.0.1): Microstructure-to-Macro Link and Cross-Market Universality

## Scope

This addendum provides a conservative theoretical framing for the empirical results in **LRD v6**.

**Claim (narrow, falsifiable):** Persistent dependence detected in Bitcoin **volatility proxies** (e.g., |log-returns|, squared returns) can be interpreted as a macroscopic manifestation of scale-invariant market microstructure mechanisms (order-flow, fragmentation, liquidity provision).

**Non-claims:** Universal predictability of Bitcoin returns, guaranteed profitability, or a single time-invariant Hurst exponent across all regimes.

## Why Volatility (Not Returns) is the Correct Target for "Memory"

Market data often exhibits an asymmetry: **returns** are frequently close to memoryless (efficient-market hypothesis), while **volatility proxies** show persistent correlations ("volatility clustering"). LRD v6 is explicitly designed to measure this asymmetry and to audit it with uncertainty estimation and falsification tests.

This separation is critical for AI applications: models assuming i.i.d. residuals in volatility estimation are structurally mis-specified when long-range dependence is present.

## Kinetic-Theory Viewpoint (Microstructure Universality)

Econophysics provides kinetic-theory approaches in which macroscopic price dynamics emerges from many interacting micro-agents (order submissions, cancellations, trades). In this view, scale-invariant microstructure can produce persistent macroscopic scaling in volatility-like state variables.

### Core References

1. **Sato, Y., & Kanazawa, K. (2025).** Strict universality of the square-root law in price impact across stocks: A complete survey of the Tokyo stock exchange. *Physical Review Letters*, 135, 257401. DOI: [10.1103/65jz-81kv](https://doi.org/10.1103/65jz-81kv)

2. **Sato, Y., & Kanazawa, K. (2024).** Does the square-root price impact law belong to the strict universal scalings?: quantitative support by a complete survey of the Tokyo stock exchange market. arXiv:2411.13965. [https://arxiv.org/abs/2411.13965](https://arxiv.org/abs/2411.13965)

3. **Bouchaud, J.-P. (2025).** The Universal Law Behind Market Price Swings. *APS Physics* (Viewpoint). [https://physics.aps.org/articles/v18/196](https://physics.aps.org/articles/v18/196)

## How This Integrates with the LRD v6 Audit Protocol

This addendum does not change the pipeline. It strengthens interpretability:

1. **Separation audit:** returns vs volatility proxies are reported separately.
2. **Uncertainty audit:** block-bootstrap CI must support H > 0.5 (for volatility proxies).
3. **Falsification audit:** phase-randomized surrogates must not reproduce the same scaling; otherwise LRD claims are rejected.
4. **Regime audit (recommended):** rolling-window DFA for Bitcoin to prevent cherry-picking and capture time variation.

## AI Relevance (Testable)

If volatility proxies are persistent, AI forecasting and risk models assuming i.i.d. noise are structurally mis-specified. LRD v6 is designed as a reproducible benchmark to test whether AI models genuinely learn persistent risk dynamics vs overfitting trends/regimes.

This has direct implications for:
- Reinforcement learning agents operating in financial environments
- Transformer-based time series models (attention mechanisms may implicitly capture LRD)
- Risk management systems calibrating VaR/CVaR under memory assumptions
