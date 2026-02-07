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


#' Produce diagnostic QQ-plot with simultaneous confidence band
#'
#' Produces a QQ-plot of the diagnostic U_i values against
#' Unif(0,1) theoretical quantiles. Linearity indicates that
#' the left-tail uniformity assumption holds. Optionally adds
#' a simultaneous 95% confidence band based on Beta order statistics.
#'
#' @param U Numeric vector of diagnostic U_i values
#' @param band Logical; add simultaneous confidence band? (default TRUE)
#' @param conf Confidence level for band (default 0.95)
#' @param B Number of Monte Carlo simulations for band calibration (default 10000)
#' @param seed Random seed for reproducibility (default 1)
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}
#' @return Invisibly, a list with components:
#'   \describe{
#'     \item{theoretical}{Theoretical quantiles}
#'     \item{observed}{Observed (sorted) U values}
#'     \item{lower}{Lower band (if band=TRUE)}
#'     \item{upper}{Upper band (if band=TRUE)}
#'     \item{outside}{Indices of observations outside band (if band=TRUE)}
#'   }
#' @export
diagnostic_qqplot <- function(U, band = TRUE, conf = 0.95, B = 10000, seed = 1, ...) {
  n <- length(U)
  if (n == 0) {
    message("No candidate T1Es detected; cannot produce QQ-plot.")
    return(invisible(NULL))
  }

  theoretical <- stats::ppoints(n)
  observed <- sort(U)

  # Compute simultaneous band if requested
  lower <- upper <- outside <- NULL
  if (band && n >= 2) {
    # Monte Carlo calibration to find pointwise level giving simultaneous coverage
    set.seed(seed)
    U_sim <- matrix(stats::runif(n * B), nrow = n, ncol = B)
    U_sim <- apply(U_sim, 2, sort)
    i <- seq_len(n)

    # Function to estimate simultaneous coverage for a given pointwise level p
    coverage_hat <- function(p) {
      tail <- (1 - p) / 2
      L <- stats::qbeta(tail, i, n + 1 - i)
      U <- stats::qbeta(1 - tail, i, n + 1 - i)
      inside <- colSums(U_sim >= L & U_sim <= U) == n
      mean(inside)
    }

    # Find p* such that simultaneous coverage = conf
    f <- function(p) coverage_hat(p) - conf
    if (f(conf) >= 0) {
      p_star <- conf
    } else {
      p_star <- stats::uniroot(f, lower = conf, upper = 0.9999, tol = 1e-4)$root
    }

    # Construct band
    tail_star <- (1 - p_star) / 2
    lower <- stats::qbeta(tail_star, i, n + 1 - i)
    upper <- stats::qbeta(1 - tail_star, i, n + 1 - i)

    # Find points outside band
    outside_ord <- which(observed < lower | observed > upper)
    outside <- order(U)[outside_ord]
  }

  # Plot
  graphics::plot(theoretical, observed,
       xlab = "Theoretical Unif(0,1) Quantiles",
       ylab = "Observed U_i Quantiles",
       main = "Diagnostic QQ-Plot",
       pch = 16, cex = 0.6, ...)
  graphics::abline(0, 1, col = "red", lwd = 2)

  if (band && n >= 2) {
    graphics::lines(theoretical, lower, col = "grey50", lty = 2)
    graphics::lines(theoretical, upper, col = "grey50", lty = 2)
    # Highlight outside points
    if (length(outside_ord) > 0) {
      graphics::points(theoretical[outside_ord], observed[outside_ord],
                       pch = 16, cex = 0.6, col = "red")
    }
  }

  invisible(list(
    theoretical = theoretical,
    observed = observed,
    lower = lower,
    upper = upper,
    outside = outside
  ))
}


#' Calibration plot (Figure 3a style)
#'
#' Multi-panel calibration plot showing observed proportion of contradictions
#' (p <= alpha AND eJAB01 > C) against alpha for a grid of C values. If the
#' method is well-calibrated, the curve should follow the reference line
#' with slope 1/up. Use \code{grid_range = c(C, C)} for a single-panel plot
#' at a specific C value.
#'
#' @param p Numeric vector of p-values (should span the full range used in
#'   estimation, not just significant results)
#' @param ejab Numeric vector of eJAB01 values (same length as \code{p})
#' @param up Upper bound for the alpha grid (default 0.1)
#' @param grid_range Length-2 numeric vector \code{c(lower, upper)} for the
#'   C grid (default \code{c(0, 1)}). Use \code{c(C, C)} for a single panel.
#' @param grid_n Number of C grid points (default 18). Panels are arranged
#'   in a 3x3 grid; more than 9 panels will span multiple pages.
#' @param n_alpha Number of alpha grid points per panel (default 200)
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}
#' @return Invisibly, a list with components \code{C_grid} (numeric vector of
#'   C values), \code{alpha_grid} (numeric vector of alpha values), and
#'   \code{proportions} (matrix with one row per C value and one column per
#'   alpha value).
#' @export
calibration_plot <- function(p, ejab, up = 0.1,
                              grid_range = c(0, 1), grid_n = 18,
                              n_alpha = 200, ...) {
  # Build C grid (same pattern as estimate_Cstar)
  if (grid_range[1] == grid_range[2]) {
    C_grid <- grid_range[1]
  } else {
    C_grid <- seq(grid_range[1], grid_range[2], length.out = grid_n)
  }

  alpha_grid <- seq(0, up, length.out = n_alpha)
  N <- sum(p < up)

  # Compute proportions: one row per C, one column per alpha
  proportions <- vapply(C_grid, function(C) {
    vapply(alpha_grid, function(a) {
      sum(p <= a & ejab > C) / N
    }, numeric(1))
  }, numeric(n_alpha))
  # vapply returns n_alpha x length(C_grid); transpose so rows = C values
  proportions <- if (length(C_grid) == 1) {
    matrix(proportions, nrow = 1)
  } else {
    t(proportions)
  }

  # Determine panel layout (max 9 per page, auto-paginates in PDF)
  n_panels <- length(C_grid)
  per_page <- min(n_panels, 9)
  n_col <- ceiling(sqrt(per_page))
  n_row <- ceiling(per_page / n_col)

  old_par <- graphics::par(mfrow = c(n_row, n_col), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old_par))

  for (i in seq_along(C_grid)) {
    prop <- proportions[i, ]
    graphics::plot(alpha_grid, prop,
         type = "l", lwd = 2,
         xlab = expression(alpha),
         ylab = "Observed Prop.",
         xlim = c(0, up),
         ylim = c(0, max(prop, up)),
         main = bquote(C == .(round(C_grid[i], 3))),
         ...)
    graphics::abline(0, 1 / up, col = "red", lty = 2, lwd = 2)
  }

  invisible(list(C_grid = C_grid, alpha_grid = alpha_grid,
                 proportions = proportions))
}


#' Full eJAB Type I Error Detection Pipeline
#'
#' Runs the complete analysis: computes eJAB01 values, estimates C*,
#' detects candidate Type I errors, and produces diagnostics.
#'
#' @param df Data frame with columns: \code{p} (p-values), \code{n}
#'   (sample sizes, > 1), \code{q} (test dimensions). An optional \code{ID}
#'   column is preserved in output.
#' @param up Upper p-value bound (default 0.05). Only results with
#'   p < up are used for estimating C*.
#' @param alpha Significance level for T1E detection (default 0.05).
#'   Must satisfy alpha <= up.
#' @param grid_range Length-2 numeric vector \code{c(lower, upper)} for the
#'   C* grid (default \code{c(1/3, 3)}). Use \code{c(1, 1)} to fix C* = 1.
#' @param grid_n Number of grid points (default 200)
#' @param plot Logical; produce calibration plot and diagnostic QQ-plot?
#'   (default TRUE)
#' @return A list with components:
#'   \describe{
#'     \item{Cstar}{Estimated threshold C*.}
#'     \item{objective}{Minimized objective function value.}
#'     \item{candidates}{Data frame of candidate T1Es with columns from input
#'       plus \code{ejab}. NULL if none detected.}
#'     \item{U}{Numeric vector of diagnostic U_i values for candidates.}
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

  # Calibration plot
  if (plot) {
    calibration_plot(df$p, ejab_vals, up = up)
  }

  # Compute diagnostics for candidates
  U <- NULL
  if (length(idx) > 0) {
    U <- diagnostic_U(df$p[idx], df$n[idx], df$q[idx], alpha, Cfit$Cstar)
    if (plot) {
      diagnostic_qqplot(U)
    }
  }

  list(
    Cstar = Cfit$Cstar,
    objective = Cfit$objective,
    candidates = if (length(idx) > 0) {
      cbind(df[idx, , drop = FALSE], ejab = ejab_vals[idx])
    } else {
      NULL
    },
    U = U,
    ejab = ejab_vals
  )
}
