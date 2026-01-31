# ejabT1E

Detects candidate Type I errors in collections of NHST results by identifying Bayes/NHST contradictions using the eJAB01 approximate objective Bayes factor.

## Installation

```r
# From this directory:
install.packages(".", repos = NULL, type = "source")

# Or with devtools:
devtools::install("path/to/ejabT1E")
```

## Usage

```r
library(ejabT1E)

df <- data.frame(
  ID = 1:1000,
  p = my_pvalues,
  n = my_sample_sizes,
  q = my_test_dimensions
)

result <- ejab_pipeline(df, up = 0.05, alpha = 0.05)

result$Cstar            # Estimated threshold
result$objective        # Minimized objective value
result$candidates       # Candidate T1Es with posterior probs and Bayes factors
result$U                # Diagnostic U_i values
result$posterior_C      # Posterior distribution over C grid
```

## Functions

| Function | Description |
|----------|-------------|
| `ejab01(p, n, q)` | Compute eJAB01 Bayes factor |
| `objective_C(C, p, ejab, up)` | Closed-form objective at a given C |
| `estimate_Cstar(p, ejab, up, grid)` | Grid search for optimal C* on [0, 10] |
| `detect_type1(p, ejab, alpha, Cstar)` | Identify candidate T1E indices |
| `diagnostic_U(p, n, q, alpha, Cstar)` | Compute U_i diagnostics (should be Unif(0,1)) |
| `diagnostic_qqplot(U)` | QQ-plot of U_i against Unif(0,1) |
| `diagnostic_qqplot_bf(U, BF)` | QQ-plot coloured by Bayes factor with colour bar |
| `posterior_t1e(p, ejab, alpha, all_objectives, grid)` | Posterior P(T1E) per result |
| `bayes_factor_t1e(posterior_prob, alpha)` | Bayes factor for T1E (posterior odds / prior odds) |
| `ejab_pipeline(df, up, alpha, grid, plot)` | Full pipeline wrapper |

## Method

A result is flagged as a candidate Type I error when:
- p-value <= alpha (NHST rejects H0), **and**
- eJAB01 > C* (Bayes factor supports H0)

C* is chosen to minimize the integrated squared deviation between the empirical contradiction rate and the theoretical T1E rate alpha/u_p across all significance levels in [0, u_p]. The grid search is over [0, 10] by default.

The posterior probability of T1E treats the objective as a negative log-likelihood under a uniform prior on C, giving a continuous probability rather than a binary classification. The Bayes factor BF = O_D / O_P compares posterior odds to prior odds (alpha / (1 - alpha)), measuring how much the data shifts the evidence for each result being a Type I error.
