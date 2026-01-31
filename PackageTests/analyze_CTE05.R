# Analysis of CTE05 clinical trial data using ejabT1E
# --------------------------------------------------

# Save console output to file
sink("CTE05_results.txt", split = TRUE)

source("../ejabT1E/R/ejab_t1e.R")

# --- Load data ---
dat <- read.csv("CTE05 (1).csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(dat), "results from CTE05\n\n")

# --- Prepare input ---
# Rename columns to match package API
df <- data.frame(
  ID        = dat$analysisId,
  p         = dat$pValue,
  n         = dat$N,
  q         = dat$q,
  # Metadata for interpretation
  nctId     = dat$nctId,
  method    = dat$statisticalMethod,
  title     = dat$briefTitle,
  outcome   = dat$outcomeTitle,
  has_lt    = dat$has_less_than
)

# --- Filter ---
# Remove rows with missing p, n, or q
complete <- complete.cases(df$p, df$n, df$q)
cat("Complete cases:", sum(complete), "/", nrow(df), "\n")
df <- df[complete, ]

# n must be > 1 (required by diagnostic_U)
df <- df[df$n > 1, ]
cat("After n > 1 filter:", nrow(df), "results\n")

# p must be in (0, 1)
df <- df[df$p > 0 & df$p < 1, ]
cat("After p in (0,1) filter:", nrow(df), "results\n")

# q must be >= 1
df <- df[df$q >= 1, ]
cat("After q >= 1 filter:", nrow(df), "results\n\n")

# --- Run pipeline ---
result <- ejab_pipeline(df, up = 0.05, alpha = 0.05, plot = FALSE)

cat("=== Results ===\n")
cat("C* estimate:", result$Cstar, "\n")
cat("Objective at C*:", result$objective, "\n")
cat("Total candidates:", nrow(result$candidates), "\n\n")

# --- Candidate summary ---
if (!is.null(result$candidates)) {
  cands <- result$candidates
  cands <- cands[order(-cands$posterior_prob), ]

  cat("Top 20 candidates by posterior P(T1E):\n")
  cat(sprintf("%-8s %-14s %8s %6s %5s %8s %8s %10s\n",
              "Rank", "NCT ID", "p-value", "n", "q", "eJAB", "P(T1E)", "BF"))
  cat(paste(rep("-", 82), collapse = ""), "\n")

  top <- head(cands, 20)
  for (i in seq_len(nrow(top))) {
    bf_str <- if (is.infinite(top$bayes_factor[i])) "Inf" else sprintf("%10.2f", top$bayes_factor[i])
    cat(sprintf("%-8d %-14s %8.4f %6d %5d %8.3f %8.3f %s\n",
                i,
                top$nctId[i],
                top$p[i],
                top$n[i],
                top$q[i],
                top$ejab[i],
                top$posterior_prob[i],
                bf_str))
  }

  cat("\n--- Posterior P(T1E) distribution ---\n")
  cat("Min:", min(cands$posterior_prob), "\n")
  cat("Median:", median(cands$posterior_prob), "\n")
  cat("Mean:", mean(cands$posterior_prob), "\n")
  cat("Max:", max(cands$posterior_prob), "\n")

  cat("\n--- Bayes Factor distribution ---\n")
  bf_finite <- cands$bayes_factor[is.finite(cands$bayes_factor)]
  cat("N finite:", length(bf_finite), "  N infinite:", sum(is.infinite(cands$bayes_factor)), "\n")
  if (length(bf_finite) > 0) {
    cat("Min:", min(bf_finite), "\n")
    cat("Median:", median(bf_finite), "\n")
    cat("Mean:", mean(bf_finite), "\n")
    cat("Max:", max(bf_finite), "\n")
  }
  cat("N with BF > 1:", sum(cands$bayes_factor > 1), "out of", nrow(cands), "\n")
}

# --- Diagnostics ---
if (!is.null(result$U)) {
  cat("\n--- Diagnostic U_i summary ---\n")
  cat("N:", length(result$U), "\n")
  cat("Range:", range(result$U), "\n")
  cat("Mean:", mean(result$U), "(expect ~0.5 under null)\n")

  # Plain QQ-plot
  pdf("CTE05_diagnostic_qqplot.pdf", width = 6, height = 6)
  diagnostic_qqplot(result$U)
  dev.off()
  cat("QQ-plot saved to CTE05_diagnostic_qqplot.pdf\n")

  # QQ-plot coloured by Bayes factor
  pdf("CTE05_diagnostic_qqplot_bf.pdf", width = 8, height = 6)
  diagnostic_qqplot_bf(result$U, BF = result$candidates$bayes_factor)
  dev.off()
  cat("BF-coloured QQ-plot saved to CTE05_diagnostic_qqplot_bf.pdf\n")
}

# --- Compare to pre-computed JAB column ---
cat("\n--- Comparison with pre-computed JAB column ---\n")
# Match back to original data by ID
matched <- match(df$ID, dat$analysisId)
jab_orig <- dat$JAB[matched]
jab_ours <- result$ejab

# Check correlation
valid <- is.finite(jab_orig) & is.finite(jab_ours)
cat("Correlation (our eJAB vs data JAB):", cor(jab_ours[valid], jab_orig[valid]), "\n")
cat("Max absolute difference:", max(abs(jab_ours[valid] - jab_orig[valid])), "\n")

# --- Save full results ---
if (!is.null(result$candidates)) {
  out <- result$candidates[order(-result$candidates$posterior_prob), ]
  write.csv(out, "CTE05_candidates.csv", row.names = FALSE)
  cat("\nFull candidate list saved to CTE05_candidates.csv\n")
}

sink()
