# ejabT1E — Package and Tests

This folder contains the R package `ejabT1E` and all associated testing and analysis materials for the eJAB Type I Error detection method.

## Structure

```
Package/
├── ejabT1E/                  R package source (installable)
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   ├── R/
│   │   └── ejab_t1e.R       All exported functions
│   ├── man/                  Rd documentation files
│   ├── .Rbuildignore
│   └── README.md             Package-specific README with API docs
│
└── PackageTests/             Testing and real-data analysis
    ├── tests/
    │   └── test_ejab.R       Unit test suite (25 tests)
    ├── analyze_CTE05.R       CTE05 clinical trial analysis script
    └── PackageOverview.ipynb Package documentation notebook
```

Files excluded by `.gitignore` (large data and reproducible output):
- `CTE05 (1).csv` — raw CTE05 dataset (~14 MB, 30,790 clinical trial results)
- `CTE05_candidates.csv` — full candidate output (~6 MB)
- `CTE05_results.txt` — console output from analysis
- `CTE05_diagnostic_qqplot.pdf` / `CTE05_diagnostic_qqplot_bf.pdf` — diagnostic plots
- `ejabT1E.Rcheck/` — R CMD check output
- `*.tar.gz`, `*.zip` — build artifacts

## Installing the package

```r
# From this directory:
install.packages("ejabT1E", repos = NULL, type = "source")

# Or with devtools:
devtools::install("ejabT1E")
```

## Running the tests

The test suite covers all exported functions with 25 tests including input validation, mathematical correctness checks, and integration tests.

```bash
cd PackageTests/tests
Rscript test_ejab.R
```

Tests verify:
- `ejab01` computation and vectorization (tasks 1.1-1.4)
- `objective_C` closed-form with gap weighting
- `estimate_Cstar` grid search on [0, 10] with minimality check
- `detect_type1` logical correctness (no false negatives or positives)
- `diagnostic_U` with varying per-test sample sizes and dimensions (task 1.2)
- `posterior_t1e` probabilities in [0, 1], posterior sums to 1 (task 1.3)
- `bayes_factor_t1e` matches O_D / O_P formula, edge cases (task 1.6)
- `diagnostic_qqplot` and `diagnostic_qqplot_bf` output structure (task 1.7)
- `ejab_pipeline` integration: all output columns present, input validation

## Running the CTE05 analysis

Requires the CTE05 dataset (`CTE05 (1).csv`) in `PackageTests/`. This file is not tracked by git due to size; obtain it from the shared data source.

```bash
cd PackageTests
Rscript analyze_CTE05.R
```

Produces:
- `CTE05_results.txt` — summary with C*, candidate counts, BF distribution
- `CTE05_candidates.csv` — full candidate list with posterior probs and Bayes factors
- `CTE05_diagnostic_qqplot.pdf` — plain QQ-plot of U_i vs Unif(0,1)
- `CTE05_diagnostic_qqplot_bf.pdf` — QQ-plot with points coloured by log10(BF)

## Exported functions

| Function | Description |
|----------|-------------|
| `ejab01(p, n, q)` | Compute eJAB01 Bayes factor |
| `objective_C(C, p, ejab, up)` | Closed-form objective for C* estimation |
| `estimate_Cstar(p, ejab, up, grid)` | Grid search for optimal C* on [0, 10] |
| `detect_type1(p, ejab, alpha, Cstar)` | Identify candidate T1E indices |
| `diagnostic_U(p, n, q, alpha, Cstar)` | Compute diagnostic U_i values |
| `diagnostic_qqplot(U)` | QQ-plot of U_i against Unif(0,1) |
| `diagnostic_qqplot_bf(U, BF)` | QQ-plot coloured by Bayes factor |
| `posterior_t1e(p, ejab, alpha, obj, grid)` | Posterior P(T1E) per result |
| `bayes_factor_t1e(posterior_prob, alpha)` | BF = posterior odds / prior odds |
| `ejab_pipeline(df, up, alpha, grid, plot)` | Full pipeline wrapper |
