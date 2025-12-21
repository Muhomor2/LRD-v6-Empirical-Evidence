# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [6.0.1] - 2025-12-19

### Added
- `ADDENDUM_THEORY.md` — Microstructure-to-macro theoretical framing linking volatility memory to scale-invariant market mechanisms. Includes AI forecasting relevance.
- `ADDENDUM_SRL_UNIVERSALITY.md` — Square-root impact law (SRL) as universality anchor. Primary sources: Sato & Kanazawa, Phys. Rev. Lett. 135, 257401 (2025); arXiv:2411.13965; APS Physics Viewpoint by Bouchaud.
- `CHANGELOG.md` — This file.
- Expanded Zenodo metadata: Description and References sections.

### Clarified
- **Narrow claim scope:** Persistent dependence is evaluated primarily on volatility proxies (|log-returns|, squared returns), not directional returns.
- **Strong audit emphasis:** Uncertainty quantification (block-bootstrap CI), surrogate falsification (phase-randomization), and scale-range sensitivity checks are mandatory validation steps.

### Unchanged
- Core empirical pipeline (DFA-2, bootstrap CI, surrogate tests)
- Data domains (Bitcoin, earthquakes, HRV, genomics)
- License (CC-BY-4.0)

## [6.0.0] - 2025-12-19

### Added
- Initial release of LRD v6: Empirical Evidence for Long-Range Dependence / Fractal Memory
- Reproducible pipeline demonstrating empirical LRD across multiple domains
- Extension of LRD v5 (Monte-Carlo validation on synthetic fGn) with real-world empirical evidence
- DFA-2 scaling estimation
- Block-bootstrap confidence intervals
- Phase-randomized surrogate tests

### Domains Covered
- Bitcoin/crypto volatility
- Earthquake activity
- Physiology (HRV)
- Genomics
