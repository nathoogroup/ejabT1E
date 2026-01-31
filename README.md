# ejabT1E — Package

This folder contains the R package `ejabT1E` for the eJAB Type I Error detection method.

## Structure

```
Package/
└── ejabT1E/                  R package source (installable)
    ├── DESCRIPTION
    ├── NAMESPACE
    ├── R/
    │   └── ejab_t1e.R       All exported functions
    ├── man/                  Rd documentation files
    ├── .Rbuildignore
    └── README.md             Package-specific README with API docs
```

Files excluded by `.gitignore`:
- `ejabT1E.Rcheck/` — R CMD check output
- `*.tar.gz`, `*.zip` — build artifacts

## Installing the package

```r
# From this directory:
install.packages("ejabT1E", repos = NULL, type = "source")

# Or with devtools:
devtools::install("ejabT1E")
```

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
