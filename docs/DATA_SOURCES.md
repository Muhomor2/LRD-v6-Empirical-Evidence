# Data Sources — LRD v6 Empirical Evidence

This document provides detailed provenance for all datasets used in the empirical evidence package.

**Important**: This package does **not** redistribute raw third-party datasets. We store only:
- Derived metrics (H/α, confidence intervals, p-values)
- Plots (DFA log-log scaling curves)
- Provenance metadata (this document)

---

## 1. Finance/Crypto — Bitcoin (BTC-USD)

### Source

**Primary**: Yahoo Finance Historical Data API

- **Endpoint**: `https://query1.finance.yahoo.com/v7/finance/download/BTC-USD`
- **Format**: CSV (OHLCV)
- **Time Range**: 2016-01-01 to 2024-12-01 (recommended)
- **Frequency**: Daily

### Query Parameters

```
period1=1451606400  # Unix timestamp for start
period2=1733011200  # Unix timestamp for end
interval=1d
events=history
includeAdjustedClose=true
```

### Series Derived

1. **Log returns**: `r_t = log(P_t) - log(P_{t-1})`
   - Expected H ≈ 0.5 (weak/no memory)
   
2. **Volatility proxy**: `|r_t|` (absolute log returns)
   - Expected H > 0.5 (volatility clustering / long memory)

### Literature Anchor

- Ding Z, Granger CWJ, Engle RF (1993). A long memory property of stock market returns and a new model. *Journal of Empirical Finance*.
- Bariviera AF (2017). The inefficiency of Bitcoin revisited: A dynamic approach. *Economics Letters*.

### License Notes

Yahoo Finance data is subject to Yahoo Terms of Service. The derived metrics and analysis outputs are original work under CC BY 4.0.

---

## 2. Geophysics — Earthquakes (USGS)

### Source

**USGS Earthquake Hazards Program** — FDSN Event Web Service

- **Endpoint**: `https://earthquake.usgs.gov/fdsnws/event/1/query`
- **Format**: GeoJSON
- **Documentation**: https://earthquake.usgs.gov/fdsnws/event/1/

### Recommended Query

```
format=geojson
starttime=2010-01-01
endtime=2023-01-01
minmagnitude=4.5
orderby=time-asc
limit=20000
```

### Series Derived

**Daily event counts**: Number of M≥4.5 earthquakes per day

Aftershock clustering creates persistent temporal structure at multiple scales, making earthquake catalogs ideal for testing long-range dependence.

### Alternative Transformations

- **Inter-event times**: Time between consecutive events
- **Magnitude series**: Sequence of event magnitudes

### License

USGS earthquake data is in the **public domain**. No redistribution restrictions apply.

### Citation

```
U.S. Geological Survey. Earthquake Hazards Program. 
FDSN Event Web Service. https://earthquake.usgs.gov/fdsnws/event/1/
```

---

## 3. Physiology — Heart Rate Variability (HRV)

### Recommended Sources

**PhysioNet** — Open databases of physiological signals

#### Open Access Databases

1. **Normal Sinus Rhythm RR Interval Database (nsrdb)**
   - URL: https://physionet.org/content/nsrdb/
   - License: Open Data Commons Attribution (ODC-By)
   - Content: RR intervals from healthy subjects

2. **MIT-BIH Normal Sinus Rhythm Database (nsrdb)**
   - URL: https://physionet.org/content/nsrdb/1.0.0/
   - License: ODC-By

3. **Congestive Heart Failure RR Interval Database (chf2db)**
   - URL: https://physionet.org/content/chf2db/
   - License: ODC-By
   - Use case: Healthy vs. pathology comparison

### Data Access

```python
import wfdb

# Example: Load record from PhysioNet
record = wfdb.rdrecord('nsrdb/16265', pn_dir='nsrdb/1.0.0')
rr_intervals = record.p_signal[:, 0]  # Extract RR intervals
```

### Series Derived

**RR deviation series**: `rr - mean(rr)` (centered RR intervals)

### Literature Anchor

- Peng CK et al. (1994). Quantification of scaling exponents and crossover phenomena in nonstationary heartbeat time series. *Chaos*.
- Goldberger AL et al. (2000). PhysioBank, PhysioToolkit, and PhysioNet. *Circulation*.

### License Notes

PhysioNet databases have explicit licenses (typically ODC-By). Always check the specific database page. Include DOI in citations.

---

## 4. Genomics — DNA Sequences (NCBI)

### Source

**NCBI GenBank / RefSeq**

- **Access**: NCBI Entrez E-utilities
- **Example Accession**: NC_000001.11 (Human chromosome 1, GRCh38.p14)
- **Format**: FASTA

### Data Access

```python
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"  # Required
handle = Entrez.efetch(db="nucleotide", id="NC_000001.11", 
                       rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
sequence = str(record.seq)
```

### Series Derived

**Purine/Pyrimidine walk**:
- A, G (purine) → +1
- C, T (pyrimidine) → -1

This binary encoding reveals long-range correlations in DNA sequences.

### Literature Anchor

- Peng CK et al. (1992). Long-range correlations in nucleotide sequences. *Nature*.
- Li W, Kaneko K (1992). Long-range correlation and partial 1/f spectrum in a noncoding DNA sequence. *Europhysics Letters*.

### License

NCBI GenBank sequences are in the **public domain**. Always record accession numbers and retrieval timestamps.

### Citation

```
NCBI Resource Coordinators (2018). Database resources of the 
National Center for Biotechnology Information. Nucleic Acids Research.
```

---

## Provenance Best Practices

For each analysis run, record:

1. **Dataset identifier**: Symbol, accession, or database name
2. **Query parameters**: Time range, filters, limits
3. **Retrieval timestamp**: When data was downloaded
4. **Transformation**: How raw data was converted to analysis series
5. **License**: Data source license terms

This information is automatically saved in `results/provenance.json` by the pipeline.

---

## Adding New Datasets

When extending this package with new domains:

1. Document the source and access method in this file
2. Note the license and any redistribution restrictions
3. Implement a loader function in `run_empirical_pipeline.py`
4. Update `provenance.json` generation
5. Add appropriate attribution and citations
