# Third Party Notices — LRD v6

This document lists third-party data sources and their licensing terms.

---

## Policy

This package follows a **derived-outputs-only** policy:

- **We do NOT redistribute raw third-party datasets**
- We redistribute only:
  - Analysis code (CC BY 4.0)
  - Derived metrics (H/α estimates, CIs, p-values)
  - Plots and figures
  - Provenance metadata

If you choose to bundle raw datasets, you must:
1. Verify redistribution is permitted
2. Include the original license
3. Provide proper attribution

---

## Data Sources and Licenses

### 1. Yahoo Finance (BTC-USD)

- **Provider**: Yahoo Finance
- **Terms**: [Yahoo Terms of Service](https://legal.yahoo.com/us/en/yahoo/terms/otos/index.html)
- **Our use**: Download via public API for analysis
- **Redistribution**: We do not redistribute raw price data

### 2. USGS Earthquake Catalog

- **Provider**: U.S. Geological Survey
- **License**: **Public Domain** (U.S. Government work)
- **Terms**: https://earthquake.usgs.gov/data/comcat/data-access.php
- **Redistribution**: Permitted without restriction
- **Citation required**: Yes (see DATA_SOURCES.md)

### 3. PhysioNet Databases

- **Provider**: PhysioNet / MIT Laboratory for Computational Physiology
- **License**: Varies by database, typically **Open Data Commons Attribution (ODC-By)**
- **Terms**: https://physionet.org/about/licenses/
- **Redistribution**: Check specific database license
- **Citation required**: Yes (include DOI)

Example databases:
- nsrdb (Normal Sinus Rhythm): ODC-By 1.0
- chf2db (CHF RR Intervals): ODC-By 1.0

### 4. NCBI GenBank

- **Provider**: National Center for Biotechnology Information
- **License**: **Public Domain**
- **Terms**: https://www.ncbi.nlm.nih.gov/home/about/policies/
- **Redistribution**: Permitted
- **Citation required**: Yes (include accession number)

---

## Software Dependencies

The following open-source packages are used:

| Package | License | URL |
|---------|---------|-----|
| NumPy | BSD-3-Clause | https://numpy.org |
| Pandas | BSD-3-Clause | https://pandas.pydata.org |
| SciPy | BSD-3-Clause | https://scipy.org |
| Matplotlib | PSF-based | https://matplotlib.org |
| Requests | Apache-2.0 | https://requests.readthedocs.io |
| Biopython | BSD-3-Clause | https://biopython.org |
| WFDB | MIT | https://github.com/MIT-LCP/wfdb-python |
| tqdm | MIT + MPL-2.0 | https://tqdm.github.io |
| statsmodels | BSD-3-Clause | https://www.statsmodels.org |

---

## Attribution Requirements

When using this package:

1. **Cite this package** (see CITATION.cff)
2. **Cite the companion package**: LRD v5 (DOI: 10.5281/zenodo.17800770)
3. **Cite original data sources** as documented in DATA_SOURCES.md

---

## Questions

For licensing questions, contact the original data providers.

For questions about this package, see the repository issues page.
