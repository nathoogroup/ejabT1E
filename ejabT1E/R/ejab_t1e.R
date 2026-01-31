#' Compute eJAB01
#'
#' Computes the approximate objective Bayes factor eJAB01 for each NHST result.
#'
#' @param p Numeric vector of p-values (0 < p < 1)
#' @param n Numeric vector of sample sizes (per test, must be > 1)
#' @param q Numeric vector of parameter dimensions (per test, >= 1)
#' @return Numeric vector of eJAB01 values
#' @examples
#' ejab01(0.03, 50, 1)
#' ejab01(c(0.01, 0.04), c(100, 200), c(1, 2))
#' @export
ejab01 <- function(p, n, q) {
  if (!all(p > 0 & p < 1)) stop("All p-values must be in (0, 1).")
  if (!all(n > 1)) stop("All sample sizes n must be > 1.")
  if (!all(q >= 1)) stop("All dimensions q must be >= 1.")
  sqrt(n) * exp(-0.5 * (n^(1/q) - 1) / n^(1/q) * stats::qchisq(1 - p, df = q))
}


#' Compute multiplicities xi_j^C
#'
#' For each unique p-value p_(j), count the number of results
#' with p-value = p_(j) AND eJAB01 > C.
#'
#' @param p_sub Numeric vector of p-values (already filtered to <= up)
#' @param ejab_sub Numeric vector of eJAB values (corresponding to p_sub)
#' @param p_unique Sorted vector of unique p-values
#' @param C Threshold value
#' @return Numeric vector of multiplicities, one per unique p-value
#' @keywords internal
compute_xi <- function(p_sub, ejab_sub, p_unique, C) {
  vapply(p_unique, function(pj) {
    sum(p_sub == pj & ejab_sub > C)
  }, numeric(1))
}


#' Closed-form objective function for C
#'
#' Evaluates the integrated squared difference between the empirical
#' contradiction rate and the theoretical T1E rate alpha/u_p.
#' Uses the closed-form expression derived from the step-function
#' structure of P_hat (see Details).
#'
#' @details
#' The objective is:
#' \deqn{\sum_{j=1}^{J} (p_{(j+1)} - p_{(j)}) \hat{P}(p_{(j)}, C)^2
#'   - \frac{1}{u_p N} \sum_{j=1}^{J} \xi_j^C (u_p^2 - p_{(j)}^2) + \frac{u_p}{3}}
#' where \eqn{p_{(J+1)} = u_p}, \eqn{\xi_j^C} counts results at p-value \eqn{p_{(j)}}
#' with eJAB > C, and \eqn{\hat{P}(p_{(j)}, C) = \frac{1}{N}\sum_{k=1}^{j}\xi_k^C}.
#'
#' @param C Candidate threshold value
#' @param p Numeric vector of all p-values
#' @param ejab Numeric vector of all eJAB01 values
#' @param up Upper bound for p-value filtering
#' @return Scalar objective function value
#' @examples
#' p <- runif(100, 0, 0.05)
#' e <- ejab01(p, rep(50, 100), rep(1, 100))
#' objective_C(1.0, p, e, up = 0.05)
#' @export
objective_C <- function(C, p, ejab, up) {
  # Filter to p < up (strict; exclude p == up to avoid boundary clumping)
  idx <- p < up
  p_sub <- p[idx]
  ejab_sub <- ejab[idx]
  N <- length(p_sub)

  if (N == 0) return(Inf)

  # Unique sorted p-values
  p_unique <- sort(unique(p_sub))
  J <- length(p_unique)

  # Compute multiplicities
  xi_C <- compute_xi(p_sub, ejab_sub, p_unique, C)

  # Cumulative P_hat at each unique p-value
  # P_hat(p_(j), C) = (1/N) * sum_{k=1}^{j} xi_k^C
  Phat_vals <- cumsum(xi_C) / N

  # Gaps: p_(j+1) - p_(j), with p_(J+1) = up
  gaps <- diff(c(p_unique, up))

  # Term 1: integral of P_hat^2
  term1 <- sum(gaps * Phat_vals^2)

  # Term 2: cross term  -2 * integral of (alpha/up) * P_hat
  term2 <- -sum(xi_C * (up^2 - p_unique^2)) / (up * N)

  # Term 3: integral of (alpha/up)^2
  term3 <- up / 3

  term1 + term2 + term3
}


#' Grid search for C*
#'
#' Finds the value of C on a grid that minimizes the objective function.
#' The default grid covers \code{[1/3, 3]} with 200 points. Pass
#' \code{grid_range = c(1, 1)} to fix C* = 1 (no search).
#'
#' @param p Numeric vector of p-values
#' @param ejab Numeric vector of eJAB01 values
#' @param up Upper bound (default 0.05)
#' @param grid_range Length-2 numeric vector \code{c(lower, upper)} specifying
#'   the grid bounds (default \code{c(1/3, 3)}). If both values are equal,
#'   C* is fixed at that value with no search.
#' @param grid_n Number of grid points (default 200). Ignored when
#'   \code{grid_range} specifies a single point.
#' @return A list with components:
#'   \describe{
#'     \item{Cstar}{The C value minimizing the objective.}
#'     \item{objective}{The minimized objective function value.}
#'     \item{all_objectives}{Numeric vector of objective values over the full grid.}
#'     \item{grid}{The grid of C values used.}
#'   }
#' @examples
#' set.seed(1)
#' p <- runif(200, 0, 0.05)
#' e <- ejab01(p, rep(50, 200), rep(1, 200))
#' fit <- estimate_Cstar(p, e, up = 0.05)
#' fit$Cstar
#' # Fix C* = 1:
#' fit1 <- estimate_Cstar(p, e, up = 0.05, grid_range = c(1, 1))
#' @export
estimate_Cstar <- function(p, ejab, up = 0.05,
                            grid_range = c(1/3, 3), grid_n = 200) {
  if (grid_range[1] == grid_range[2]) {
    grid <- grid_range[1]
  } else {
    grid <- seq(grid_range[1], grid_range[2], length.out = grid_n)
  }
  vals <- vapply(grid, objective_C, numeric(1), p = p, ejab = ejab, up = up)

  idx <- which.min(vals)
  list(
    Cstar = grid[idx],
    objective = vals[idx],
    all_objectives = vals,
    grid = grid
  )
}


#' Detect candidate Type I errors
#'
#' Identifies results that are Bayes/NHST contradictions: p-value <= alpha
#' (NHST rejects) and eJAB01 > Cstar (Bayes factor supports H0).
#'
#' @param p Numeric vector of p-values
#' @param ejab Numeric vector of eJAB01 values
#' @param alpha Significance level (must be <= up used in estimation)
#' @param Cstar Estimated threshold
#' @return Integer vector of indices of candidate T1Es
#' @export
detect_type1 <- function(p, ejab, alpha, Cstar) {
  which(p < alpha & ejab > Cstar)
}


#' Compute diagnostic U_i values
#'
#' For each candidate T1E, computes a diagnostic U_i that should follow
#' Unif(0,1) if the left-tail uniformity assumption holds and the
#' selected results are true Type I errors.
#'
#' @details
#' The diagnostic is defined as:
#' \deqn{U_i = \frac{p_i - d_i}{\alpha - d_i}}
#' where
#' \deqn{d_i = 1 - F_{\chi^2_{q_i}}\left(\frac{2 n_i^{1/q_i}}{n_i^{1/q_i} - 1}
#'   \log\frac{\sqrt{n_i}}{C^*}\right)}
#'
#' @param p Numeric vector of p-values (for detected T1Es only)
#' @param n Numeric vector of sample sizes (per test, for detected T1Es only; must be > 1)
#' @param q Numeric vector of dimensions (per test, for detected T1Es only)
#' @param alpha Significance level
#' @param Cstar Estimated threshold
#' @return Numeric vector of U_i values
#' @export
diagnostic_U <- function(p, n, q, alpha, Cstar) {
  if (any(n <= 1)) stop("Sample sizes n must be > 1 for diagnostic computation.")
  d <- 1 - stats::pchisq(
    (2 * n^(1/q) / (n^(1/q) - 1)) * log(sqrt(n) / Cstar),
    df = q
  )
  (p - d) / (alpha - d)
}


#' Produce diagnostic QQ-plot
#'
#' Produces a QQ-plot of the diagnostic U_i values against
#' Unif(0,1) theoretical quantiles. Linearity indicates that
#' the left-tail uniformity assumption holds.
#'
#' @param U Numeric vector of diagnostic U_i values
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}
#' @return Invisibly, a list with components \code{theoretical} and \code{observed}.
#' @export
diagnostic_qqplot <- function(U, ...) {
  T_len <- length(U)
  if (T_len == 0) {
    message("No candidate T1Es detected; cannot produce QQ-plot.")
    return(invisible(NULL))
  }
  theoretical <- stats::qunif(stats::ppoints(T_len))
  observed <- sort(U)

  graphics::plot(theoretical, observed,
       xlab = "Theoretical Unif(0,1) Quantiles",
       ylab = "Observed U_i Quantiles",
       main = "Diagnostic QQ-Plot",
       pch = 16, cex = 0.6, ...)
  graphics::abline(0, 1, col = "red", lwd = 2)
  invisible(list(theoretical = theoretical, observed = observed))
}


#' Compute posterior probability of Type I error for each result
#'
#' Treats C* as the posterior mode under a uniform prior on the grid.
#' Converts objective function values to a posterior distribution over C,
#' then for each result computes P(T1E) by summing posterior mass over
#' all C values satisfying the contradiction criterion.
#'
#' @details
#' The posterior over C is obtained by:
#' \deqn{\pi(C_k | \text{data}) \propto \exp\{-\text{objective}(C_k)\}}
#' assuming a uniform prior over the grid. For each result i, the
#' posterior probability of being a Type I error is:
#' \deqn{P(\text{T1E}_i) = \sum_{k: p_i \le \alpha,\; \text{eJAB}_i > C_k} \pi(C_k | \text{data})}
#'
#' @param p Numeric vector of p-values (for results of interest)
#' @param ejab Numeric vector of eJAB01 values
#' @param alpha Significance level
#' @param all_objectives Numeric vector of objective function values over the grid
#' @param grid Numeric vector of C values corresponding to all_objectives
#' @return Numeric vector of posterior T1E probabilities, one per result
#' @export
posterior_t1e <- function(p, ejab, alpha, all_objectives, grid) {
  # Convert objectives to unnormalized posterior: exp(-obj)
  # Subtract max of log-values for numerical stability
  log_unnorm <- -all_objectives
  log_unnorm <- log_unnorm - max(log_unnorm)
  unnorm <- exp(log_unnorm)
  posterior_C <- unnorm / sum(unnorm)

  # For each result, sum posterior(c_k) over all c_k where
  # p_i <= alpha AND eJAB_i > c_k
  vapply(seq_along(p), function(i) {
    if (p[i] >= alpha) return(0)
    qualifying <- ejab[i] > grid
    sum(posterior_C[qualifying])
  }, numeric(1))
}


#' Compute Bayes factor for Type I error
#'
#' For each candidate T1E, computes a Bayes factor measuring the evidence
#' that the result is a Type I error, relative to the prior.
#'
#' @details
#' The posterior odds are:
#' \deqn{O_D = \frac{P(\text{T1E} \mid \text{Data})}{1 - P(\text{T1E} \mid \text{Data})}}
#' The prior odds (using alpha as the prior probability of a T1E) are:
#' \deqn{O_P = \frac{\alpha}{1 - \alpha}}
#' The Bayes factor is \eqn{BF = O_D / O_P}.
#'
#' @param posterior_prob Numeric vector of posterior T1E probabilities
#'   (from \code{\link{posterior_t1e}})
#' @param alpha Significance level (prior probability of T1E)
#' @return Numeric vector of Bayes factors. Values > 1 indicate data increased
#'   the evidence for T1E relative to the prior. Inf if posterior_prob = 1.
#' @export
bayes_factor_t1e <- function(posterior_prob, alpha) {
  if (any(posterior_prob < 0 | posterior_prob > 1)) {
    stop("posterior_prob must be in [0, 1].")
  }
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1).")
  prior_odds <- alpha / (1 - alpha)
  posterior_odds <- posterior_prob / (1 - posterior_prob)
  posterior_odds / prior_odds
}


#' Produce diagnostic QQ-plot coloured by Bayes factor
#'
#' Produces a QQ-plot of the diagnostic U_i values against
#' Unif(0,1) theoretical quantiles, with points coloured by
#' their Bayes factor for Type I error. Includes a colour bar.
#'
#' @param U Numeric vector of diagnostic U_i values
#' @param BF Numeric vector of Bayes factors (same length as U).
#'   If NULL, all points are plotted in black (no colour mapping).
#' @param log_BF Logical; if TRUE (default), colour scale uses log10(BF)
#' @param col_low Colour for lowest BF values (default "blue")
#' @param col_high Colour for highest BF values (default "red")
#' @param n_cols Number of colours in the gradient (default 256)
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}
#' @return Invisibly, a list with components \code{theoretical}, \code{observed},
#'   and \code{BF}.
#' @export
diagnostic_qqplot_bf <- function(U, BF = NULL, log_BF = TRUE,
                                 col_low = "blue", col_high = "red",
                                 n_cols = 256, ...) {
  T_len <- length(U)
  if (T_len == 0) {
    message("No candidate T1Es detected; cannot produce QQ-plot.")
    return(invisible(NULL))
  }
  theoretical <- stats::qunif(stats::ppoints(T_len))
  ord <- order(U)
  observed <- U[ord]

  if (!is.null(BF)) {
    BF_sorted <- BF[ord]
    if (log_BF) {
      val <- log10(pmax(BF_sorted, .Machine$double.eps))
    } else {
      val <- BF_sorted
    }
    rng <- range(val, finite = TRUE)
    if (rng[1] == rng[2]) rng <- rng + c(-0.5, 0.5)
    scaled <- (val - rng[1]) / (rng[2] - rng[1])
    scaled[scaled < 0] <- 0
    scaled[scaled > 1] <- 1

    pal <- grDevices::colorRampPalette(c(col_low, col_high))(n_cols)
    cols <- pal[pmax(1, ceiling(scaled * n_cols))]
  } else {
    cols <- "black"
    BF_sorted <- NULL
  }

  # Save and restore par
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  if (!is.null(BF)) {
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
    graphics::par(mar = c(5, 4, 4, 1))
  }

  graphics::plot(theoretical, observed,
       xlab = "Theoretical Unif(0,1) Quantiles",
       ylab = "Observed U_i Quantiles",
       main = "Diagnostic QQ-Plot",
       col = cols, pch = 16, cex = 0.6, ...)
  graphics::abline(0, 1, col = "grey40", lwd = 2)

  # Colour bar
  if (!is.null(BF)) {
    graphics::par(mar = c(5, 0.5, 4, 3))
    bar_vals <- seq(0, 1, length.out = n_cols)
    pal <- grDevices::colorRampPalette(c(col_low, col_high))(n_cols)
    graphics::image(1, bar_vals, matrix(bar_vals, nrow = 1), col = pal,
          axes = FALSE, xlab = "", ylab = "")
    labels_at <- seq(0, 1, length.out = 5)
    labels_val <- round(rng[1] + labels_at * (rng[2] - rng[1]), 2)
    graphics::axis(4, at = labels_at, labels = labels_val, las = 1, cex.axis = 0.8)
    graphics::mtext(if (log_BF) expression(log[10](BF)) else "BF",
          side = 4, line = 2, cex = 0.8)
  }

  invisible(list(theoretical = theoretical, observed = observed, BF = BF_sorted))
}


#' Full eJAB Type I Error Detection Pipeline
#'
#' Runs the complete analysis: computes eJAB01 values, estimates C*,
#' detects candidate Type I errors, produces diagnostics, and computes
#' posterior probabilities.
#'
#' @param df Data frame with columns: \code{p} (p-values), \code{n}
#'   (sample sizes, > 1), \code{q} (test dimensions). An optional \code{ID}
#'   column is preserved in output.
#' @param up Upper p-value bound (default 0.05). Only results with
#'   p <= up are used for estimating C*.
#' @param alpha Significance level for T1E detection (default 0.05).
#'   Must satisfy alpha <= up.
#' @param grid_range Length-2 numeric vector \code{c(lower, upper)} for the
#'   C* grid (default \code{c(1/3, 3)}). Use \code{c(1, 1)} to fix C* = 1.
#' @param grid_n Number of grid points (default 200)
#' @param plot Logical; produce diagnostic QQ-plot? (default TRUE)
#' @return A list with components:
#'   \describe{
#'     \item{Cstar}{Estimated threshold C*.}
#'     \item{objective}{Minimized objective function value.}
#'     \item{candidates}{Data frame of candidate T1Es with columns from input
#'       plus \code{ejab}, \code{posterior_prob}, and \code{bayes_factor}.
#'       NULL if none detected.}
#'     \item{U}{Numeric vector of diagnostic U_i values for candidates.}
#'     \item{posterior_C}{Data frame with columns \code{C} and \code{posterior}
#'       giving the posterior distribution over the grid.}
#'     \item{ejab}{Numeric vector of eJAB01 values for all input results.}
#'   }
#' @examples
#' set.seed(42)
#' df <- data.frame(
#'   ID = 1:100,
#'   p = runif(100, 0, 0.05),
#'   n = sample(20:200, 100, replace = TRUE),
#'   q = rep(1, 100)
#' )
#' result <- ejab_pipeline(df, up = 0.05, alpha = 0.05, plot = FALSE)
#' result$Cstar
#' head(result$candidates)
#' @export
ejab_pipeline <- function(df, up = 0.05, alpha = 0.05,
                          grid_range = c(1/3, 3), grid_n = 200,
                          plot = TRUE) {
  # Input validation
  if (!all(c("p", "n", "q") %in% names(df))) {
    stop("df must contain columns 'p', 'n', and 'q'.")
  }
  if (alpha > up) {
    stop("alpha must be <= up. The method requires alpha in [0, up].")
  }
  if (any(df$n <= 1)) {
    stop("All sample sizes (n) must be > 1.")
  }

  # Compute eJAB01
  ejab_vals <- ejab01(df$p, df$n, df$q)

  # Estimate C*
  Cfit <- estimate_Cstar(df$p, ejab_vals, up, grid_range, grid_n)

  # Detect candidate T1Es
  idx <- detect_type1(df$p, ejab_vals, alpha, Cfit$Cstar)

  # Compute diagnostics for candidates
  U <- NULL
  post_prob <- NULL
  bf <- NULL

  if (length(idx) > 0) {
    U <- diagnostic_U(df$p[idx], df$n[idx], df$q[idx], alpha, Cfit$Cstar)

    # Posterior probabilities
    post_prob <- posterior_t1e(
      df$p[idx], ejab_vals[idx], alpha,
      Cfit$all_objectives, Cfit$grid
    )

    # Bayes factors
    bf <- bayes_factor_t1e(post_prob, alpha)

    if (plot) {
      diagnostic_qqplot_bf(U, BF = bf)
    }
  }

  # Posterior distribution over C
  log_unnorm <- -Cfit$all_objectives
  log_unnorm <- log_unnorm - max(log_unnorm)
  unnorm <- exp(log_unnorm)
  posterior_C <- unnorm / sum(unnorm)

  list(
    Cstar = Cfit$Cstar,
    objective = Cfit$objective,
    candidates = if (length(idx) > 0) {
      cbind(df[idx, , drop = FALSE],
            ejab = ejab_vals[idx],
            posterior_prob = post_prob,
            bayes_factor = bf)
    } else {
      NULL
    },
    U = U,
    posterior_C = data.frame(C = Cfit$grid, posterior = posterior_C),
    ejab = ejab_vals
  )
}
