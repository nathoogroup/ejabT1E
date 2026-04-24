# ejabT1E

Detects candidate Type I errors in collections of NHST results by identifying Bayes/NHST contradictions using the eJAB01 approximate objective Bayes factor.

## Installation

```r
# From this directory:
install.packages(".", repos = NULL, type = "source")

# Or with devtools:
devtools::install("path/to/ejabT1E")
```

## Quick Start

```r
library(ejabT1E)

# Your data needs: p-values, sample sizes, parameter dimensions
df <- data.frame(
  p = my_pvalues,   # p-values in (0, 1)
  n = my_sample_sizes,  # sample sizes > 1
  q = my_test_dimensions  # parameter dimensions >= 1 (usually 1 for t-tests)
)

# Run full pipeline
result <- ejab_pipeline(df, up = 0.05, alpha = 0.05)

# Results
result$Cstar       # Estimated threshold C*
result$objective   # Minimized objective value
result$candidates  # Data frame of candidate T1Es
result$U           # Diagnostic U values (should be Unif(0,1))
result$ejab        # eJAB01 values for all results
```

## Functions

| Function | Description |
|----------|-------------|
| `ejab01(p, n, q)` | Compute eJAB01 Bayes factor |
| `objective_C(C, p, ejab, up)` | Evaluate objective function at a given C |
| `estimate_Cstar(p, ejab, up, grid_range, grid_n)` | Grid search for optimal C* |
| `detect_type1(p, ejab, alpha, Cstar)` | Identify candidate T1E indices |
| `diagnostic_U(p, n, q, alpha, Cstar)` | Compute diagnostic U values |
| `diagnostic_qqplot(U, alpha, Cstar, band, conf, B, seed)` | QQ-plot with best-fit line and simultaneous confidence band |
| `calibration_plot(p, ejab, up, alpha, grid_range, grid_n, n_alpha, n, q)` | Adaptive calibration plot using C*(alpha) (3 plots) |
| `ejab_pipeline(df, up, alpha, grid_range, grid_n, plot)` | Full analysis pipeline |

## Method

A result is flagged as a candidate Type I error when:
- p-value < alpha (NHST rejects H0), **and**
- eJAB01 > C* (Bayes factor supports H0)

### C* Estimation

C* is chosen to minimize the integrated squared deviation between the empirical contradiction rate and the theoretical T1E rate. Default grid: [1/3, 3] with 200 points.

```r
# Use default grid
fit <- estimate_Cstar(p, ejab, up = 0.05)

# Fix C* = 1 (no estimation)
fit <- estimate_Cstar(p, ejab, up = 0.05, grid_range = c(1, 1))

# Custom grid
fit <- estimate_Cstar(p, ejab, up = 0.05, grid_range = c(0.5, 2), grid_n = 100)
```

### Diagnostic QQ-Plot

The diagnostic U values should follow Unif(0,1) if detected candidates are true Type I errors. The QQ-plot overlays a best-fit line (OLS) instead of the 45-degree line, since C* is estimated. Includes a 95% simultaneous confidence band (Monte Carlo calibrated), transformed to follow the best-fit line.

```r
# With confidence band and title showing alpha/C*
diagnostic_qqplot(U, alpha = 0.05, Cstar = 0.01, band = TRUE, conf = 0.95)

# Without band
diagnostic_qqplot(U, band = FALSE)
```

Points outside the band are highlighted in red and returned in `$outside`.

### Calibration Plot

Adaptive calibration plot using C*(alpha). For each alpha on a grid in [0, up], finds C*(alpha) such that the observed proportion of contradictions is closest to alpha/up. Produces three plots:

1. **Calibration curve** using the adaptive C*(alpha) — should track the reference line (slope 1/up) closely
2. **C*(alpha) vs alpha** — shows how the threshold varies (expected: decreasing)
3. **Diagnostic QQ-plot** at the specified alpha using C*(alpha)

```r
# Default (3 plots, QQ-plot at alpha = 0.05)
calibration_plot(p, ejab, up = 0.1, alpha = 0.05, n = n_vec, q = q_vec)

# Without QQ-plot (omit n and q)
calibration_plot(p, ejab, up = 0.1)
```

## Input Requirements

- **p**: p-values strictly in (0, 1)
- **n**: sample sizes > 1
- **q**: parameter dimensions >= 1
- **up**: upper bound for p-values used in C* estimation (default 0.05)
- **alpha**: significance level for T1E detection, must be <= up (default 0.05)

## Example

```r
library(ejabT1E)

# Simulated data
set.seed(42)
df <- data.frame(
  p = runif(500, 0, 0.05),
  n = sample(50:200, 500, replace = TRUE),
  q = rep(1, 500)
)

# Run analysis
result <- ejab_pipeline(df, up = 0.05, alpha = 0.05, plot = TRUE)

cat("C* =", result$Cstar, "\n")
cat("Candidates:", nrow(result$candidates), "\n")
cat("Outside band:", length(result$U[result$U < 0 | result$U > 1]), "\n")
```

## Reference

Velidi P, Wei Z, Kalaria SN, Liu Y, Laumont CM, Nelson BH, Nathoo FS. (2025). Generalized Jeffreys's approximate objective Bayes factor: Model-selection consistency, finite-sample accuracy and statistical evidence in 71,126 clinical trial findings. PNAS.
