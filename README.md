# Chapter3-surrogacyTLS

Repository of analysis code for Chapter 3 of my PhD thesis: **"RISE: rank-based identification of high-dimensional surrogate markers in the single-trial setting"**.

RISE is applied here to transcriptomic data from vaccine trials, to identify high-dimensional gene expression (and gene set) surrogate markers for clinical/immunological outcomes, using the [`SurrogateRank`](https://github.com/laylaparast/SurrogateRank) package's `rise.screen` / `rise.evaluate` framework.

## Repository structure

```
analysis/
├── tutorial/         # Self-contained worked examples of the RISE method on simulated data
├── preprocessing/     # Scripts to clean, harmonise, and merge the raw trial datasets
├── descriptive/       # Descriptive/exploratory analyses of the harmonised datasets
└── application/       # Application of RISE to identify surrogate markers in real trial data
```

### `analysis/tutorial`
Illustrative, standalone scripts demonstrating the RISE workflow (screening + evaluation, sample-splitting) on simulated data, plus a bivariate normal illustration of the underlying theory.

### `analysis/preprocessing`
Scripts for reading in and cleaning the individual study datasets (Hamburg, EBOVAC2, PREVAC, IS2 clinical/expression/immune-response data), harmonising clinical variables across studies, merging datasets, and preparing gene set definitions (BTM, BG3M). `preprocessing_master.R` sources these in the required order to reproduce the full pipeline from raw data to an analysis-ready, harmonised dataset.

### `analysis/descriptive`
Exploratory/descriptive analyses of the preprocessed data (e.g. Ebolavirus and SDY1276 datasets).

### `analysis/application`
Application of the RISE method to the harmonised trial data, including study-specific scripts (RISE-Ad26MVA, RISE-rVSV, RISE-SDY1276) and shared RISE signature code. `application_master.R` sources these scripts to run the full applied analysis.

## Requirements

R with the following packages (non-exhaustive, see individual scripts for full dependencies):
- [`SurrogateRank`](https://github.com/laylaparast/SurrogateRank)
- `fs`

## Usage

1. Run `analysis/preprocessing/preprocessing_master.R` to reproduce the harmonised analysis dataset from raw sources.
2. Run `analysis/application/application_master.R` to reproduce the applied RISE surrogate marker analysis.
3. See `analysis/tutorial/` for a minimal, self-contained introduction to the RISE method before diving into the applied scripts.
