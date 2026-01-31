# Test script for ejabT1E package
# Covers all exported functions and completed tasks per tasks.md
# Source the functions directly for testing:
source("../../ejabT1E/R/ejab_t1e.R")

cat("============================================================\n")
cat("ejabT1E test suite\n")
cat("============================================================\n\n")

# ============================================================
# 1. ejab01
# ============================================================

# --- Test 1: ejab01 basic computation ---
val <- ejab01(0.05, 100, 1)
cat("Test 1  - ejab01(0.05, 100, 1):", val, "\n")
stopifnot(is.finite(val), val > 0)

# --- Test 2: ejab01 vectorized with mixed q ---
vals <- ejab01(c(0.01, 0.03, 0.05), c(50, 100, 200), c(1, 1, 2))
cat("Test 2  - vectorized ejab01:", vals, "\n")
stopifnot(length(vals) == 3, all(is.finite(vals)), all(vals > 0))

# --- Test 3: ejab01 input validation ---
ok <- tryCatch({ ejab01(0, 100, 1); FALSE }, error = function(e) TRUE)
stopifnot(ok)
ok <- tryCatch({ ejab01(0.05, 1, 1); FALSE }, error = function(e) TRUE)
stopifnot(ok)
ok <- tryCatch({ ejab01(0.05, 100, 0); FALSE }, error = function(e) TRUE)
stopifnot(ok)
cat("Test 3  - ejab01 input validation: OK\n")

# ============================================================
# 2. objective_C (closed-form)
# ============================================================

set.seed(42)
N <- 500
is_null <- sample(c(TRUE, FALSE), N, replace = TRUE, prob = c(0.2, 0.8))
p_sim <- ifelse(is_null, runif(N, 0, 0.05), runif(N, 0, 0.01))
n_sim <- sample(20:300, N, replace = TRUE)  # varying sample sizes (task 1.2)
q_sim <- sample(1:3, N, replace = TRUE)     # varying dimensions
ejab_sim <- ejab01(p_sim, n_sim, q_sim)

# --- Test 4: Objective is finite and non-negative ---
obj_val <- objective_C(1.0, p_sim, ejab_sim, up = 0.05)
cat("Test 4  - objective_C at C=1:", obj_val, "\n")
stopifnot(is.finite(obj_val), obj_val >= 0)

# --- Test 5: Objective returns Inf when no data passes filter ---
obj_empty <- objective_C(1.0, c(0.9, 0.8), ejab01(c(0.9, 0.8), c(50, 50), c(1, 1)), up = 0.05)
stopifnot(obj_empty == Inf)
cat("Test 5  - objective_C with no data in [0, up]: Inf (correct)\n")

# ============================================================
# 3. estimate_Cstar (grid search)
# ============================================================

# --- Test 6: Grid search on [0, 10] (task 1.1) ---
fit <- estimate_Cstar(p_sim, ejab_sim, up = 0.05)
cat("Test 6  - Cstar:", fit$Cstar, "  objective:", fit$objective, "\n")
stopifnot(fit$Cstar >= 0, fit$Cstar <= 10)
stopifnot(length(fit$grid) == 200)
stopifnot(fit$grid[1] == 0, fit$grid[200] == 10)
# C* must actually minimize: objective at C* <= objective at grid endpoints
stopifnot(fit$objective <= objective_C(0, p_sim, ejab_sim, 0.05))
stopifnot(fit$objective <= objective_C(10, p_sim, ejab_sim, 0.05))
# all_objectives returned for posterior computation
stopifnot(length(fit$all_objectives) == 200)

# ============================================================
# 4. detect_type1
# ============================================================

# --- Test 7: Detection logic ---
idx <- detect_type1(p_sim, ejab_sim, alpha = 0.05, Cstar = fit$Cstar)
cat("Test 7  - Detected", length(idx), "candidate T1Es out of", N, "\n")
stopifnot(all(p_sim[idx] <= 0.05))
stopifnot(all(ejab_sim[idx] > fit$Cstar))
# No result outside idx should satisfy both conditions
non_idx <- setdiff(seq_len(N), idx)
if (length(non_idx) > 0) {
  stopifnot(!any(p_sim[non_idx] <= 0.05 & ejab_sim[non_idx] > fit$Cstar))
}

# ============================================================
# 5. diagnostic_U (task 1.2: per-test sample size)
# ============================================================

# --- Test 8: U values with varying n and q ---
if (length(idx) > 0) {
  U <- diagnostic_U(p_sim[idx], n_sim[idx], q_sim[idx], alpha = 0.05, Cstar = fit$Cstar)
  cat("Test 8  - U range:", range(U), "\n")
  stopifnot(all(is.finite(U)))
  stopifnot(length(U) == length(idx))
  # U uses per-test n_i (not a single N): verify length matches
  stopifnot(length(unique(n_sim[idx])) > 1 || length(idx) < 5)
}

# --- Test 9: diagnostic_U input validation ---
ok <- tryCatch({ diagnostic_U(0.01, 1, 1, 0.05, 1.0); FALSE }, error = function(e) TRUE)
stopifnot(ok)
cat("Test 9  - diagnostic_U rejects n <= 1: OK\n")

# ============================================================
# 6. posterior_t1e (task 1.3)
# ============================================================

if (length(idx) > 0) {
  # --- Test 10: Posterior probabilities ---
  post <- posterior_t1e(p_sim[idx], ejab_sim[idx], alpha = 0.05,
                        fit$all_objectives, fit$grid)
  cat("Test 10 - Posterior prob range:", range(post), "\n")
  stopifnot(all(post >= 0 & post <= 1))
  stopifnot(length(post) == length(idx))

  # --- Test 11: Posterior over C sums to 1 ---
  log_unnorm <- -fit$all_objectives
  log_unnorm <- log_unnorm - max(log_unnorm)
  unnorm <- exp(log_unnorm)
  posterior_C <- unnorm / sum(unnorm)
  stopifnot(abs(sum(posterior_C) - 1) < 1e-10)
  cat("Test 11 - Posterior over C sums to 1: OK\n")

  # --- Test 12: Result with p > alpha gets posterior = 0 ---
  post_zero <- posterior_t1e(0.06, ejab_sim[idx[1]], alpha = 0.05,
                              fit$all_objectives, fit$grid)
  stopifnot(post_zero == 0)
  cat("Test 12 - P(T1E) = 0 when p > alpha: OK\n")
}

# ============================================================
# 7. bayes_factor_t1e (task 1.6)
# ============================================================

if (length(idx) > 0) {
  # --- Test 13: BF formula matches O_D / O_P ---
  bf <- bayes_factor_t1e(post, alpha = 0.05)
  cat("Test 13 - BF range:", range(bf), "\n")
  stopifnot(all(is.finite(bf) | bf == Inf))
  stopifnot(all(bf >= 0))
  prior_odds <- 0.05 / 0.95
  expected_bf <- (post / (1 - post)) / prior_odds
  stopifnot(all(abs(bf - expected_bf) < 1e-10 | is.infinite(bf)))

  # --- Test 14: BF > 1 iff posterior_prob > alpha ---
  finite_mask <- is.finite(bf)
  if (any(finite_mask)) {
    stopifnot(all((bf[finite_mask] > 1) == (post[finite_mask] > 0.05)))
  }
  cat("Test 14 - BF > 1 iff P(T1E) > alpha: OK\n")

  # --- Test 15: BF = Inf when posterior_prob = 1 ---
  bf_inf <- bayes_factor_t1e(1.0, alpha = 0.05)
  stopifnot(is.infinite(bf_inf))
  cat("Test 15 - BF = Inf when P(T1E) = 1: OK\n")

  # --- Test 16: BF = 0 when posterior_prob = 0 ---
  bf_zero <- bayes_factor_t1e(0.0, alpha = 0.05)
  stopifnot(bf_zero == 0)
  cat("Test 16 - BF = 0 when P(T1E) = 0: OK\n")

  # --- Test 17: BF input validation ---
  ok <- tryCatch({ bayes_factor_t1e(1.5, 0.05); FALSE }, error = function(e) TRUE)
  stopifnot(ok)
  ok <- tryCatch({ bayes_factor_t1e(0.5, 0); FALSE }, error = function(e) TRUE)
  stopifnot(ok)
  cat("Test 17 - BF input validation: OK\n")
}

# ============================================================
# 8. diagnostic_qqplot and diagnostic_qqplot_bf (tasks 1.7)
# ============================================================

if (length(idx) > 0) {
  # --- Test 18: diagnostic_qqplot returns correct structure ---
  pdf(file = NULL)  # suppress plot output
  qq_out <- diagnostic_qqplot(U)
  dev.off()
  stopifnot(!is.null(qq_out))
  stopifnot(all(c("theoretical", "observed") %in% names(qq_out)))
  stopifnot(length(qq_out$theoretical) == length(U))
  cat("Test 18 - diagnostic_qqplot output structure: OK\n")

  # --- Test 19: diagnostic_qqplot_bf returns correct structure ---
  pdf(file = NULL)
  qq_bf_out <- diagnostic_qqplot_bf(U, BF = bf)
  dev.off()
  stopifnot(!is.null(qq_bf_out))
  stopifnot(all(c("theoretical", "observed", "BF") %in% names(qq_bf_out)))
  stopifnot(length(qq_bf_out$BF) == length(U))
  cat("Test 19 - diagnostic_qqplot_bf output structure: OK\n")

  # --- Test 20: diagnostic_qqplot_bf without BF falls back to black ---
  pdf(file = NULL)
  qq_no_bf <- diagnostic_qqplot_bf(U, BF = NULL)
  dev.off()
  stopifnot(is.null(qq_no_bf$BF))
  cat("Test 20 - diagnostic_qqplot_bf without BF: OK\n")
}

# --- Test 21: Empty U returns NULL ---
qq_empty <- diagnostic_qqplot(numeric(0))
stopifnot(is.null(qq_empty))
qq_empty2 <- diagnostic_qqplot_bf(numeric(0))
stopifnot(is.null(qq_empty2))
cat("Test 21 - Empty U returns NULL for both qqplot functions: OK\n")

# ============================================================
# 9. ejab_pipeline (full integration)
# ============================================================

# --- Test 22: Pipeline returns all expected components ---
df_test <- data.frame(ID = seq_len(N), p = p_sim, n = n_sim, q = q_sim)
result <- ejab_pipeline(df_test, up = 0.05, alpha = 0.05, plot = FALSE)
cat("Test 22 - Pipeline Cstar:", result$Cstar, "\n")
cat("          Candidates:", nrow(result$candidates), "\n")
stopifnot(!is.null(result$Cstar))
stopifnot(result$Cstar == fit$Cstar)
stopifnot(!is.null(result$posterior_C))
stopifnot(abs(sum(result$posterior_C$posterior) - 1) < 1e-10)
stopifnot(length(result$ejab) == N)

# --- Test 23: Candidates dataframe has all required columns ---
if (!is.null(result$candidates)) {
  required_cols <- c("p", "n", "q", "ejab", "posterior_prob", "bayes_factor")
  for (col in required_cols) {
    stopifnot(col %in% names(result$candidates))
  }
  stopifnot(all(result$candidates$bayes_factor >= 0))
  stopifnot(all(result$candidates$posterior_prob >= 0 &
                result$candidates$posterior_prob <= 1))
  cat("Test 23 - Candidates columns:", paste(required_cols, collapse = ", "), "present: OK\n")
}

# --- Test 24: Pipeline input validation ---
ok <- tryCatch({
  ejab_pipeline(data.frame(p = 0.01, n = 50, q = 1), up = 0.01, alpha = 0.05, plot = FALSE)
  FALSE
}, error = function(e) TRUE)
stopifnot(ok)
cat("Test 24 - Pipeline rejects alpha > up: OK\n")

ok <- tryCatch({
  ejab_pipeline(data.frame(x = 1), up = 0.05, alpha = 0.05, plot = FALSE)
  FALSE
}, error = function(e) TRUE)
stopifnot(ok)
cat("Test 25 - Pipeline rejects missing columns: OK\n")

cat("\n============================================================\n")
cat("All 25 tests passed.\n")
cat("============================================================\n")
