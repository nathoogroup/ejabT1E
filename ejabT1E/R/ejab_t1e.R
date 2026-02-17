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
#' the left-tail uniformity assumption holds. A best-fit line
#' (via OLS) is overlaid instead of the 45-degree line, since
#' C* is estimated. Optionally adds a simultaneous confidence
#' band based on Beta order statistics, transformed to follow the
#' best-fit line rather than the 45-degree line.
#'
#' @param U Numeric vector of diagnostic U_i values
#' @param alpha Significance level used for detection (displayed in title).
#'   If NULL, omitted from title.
#' @param Cstar C* value used for detection (displayed in title).
#'   If NULL, omitted from title.
#' @param band Logical; add simultaneous confidence band? (default TRUE)
#' @param conf Confidence level for band (default 0.95)
#' @param B Number of Monte Carlo simulations for band calibration (default 10000)
#' @param seed Random seed for reproducibility (default 1)
#' @param ... Additional arguments (currently unused, kept for compatibility)
#' @return Invisibly, a list with components:
#'   \describe{
#'     \item{theoretical}{Theoretical quantiles}
#'     \item{observed}{Observed (sorted) U values}
#'     \item{intercept}{OLS intercept}
#'     \item{slope}{OLS slope}
#'     \item{lower}{Lower band (if band=TRUE)}
#'     \item{upper}{Upper band (if band=TRUE)}
#'     \item{outside}{Indices of observations outside band (if band=TRUE)}
#'   }
#' @export
diagnostic_qqplot <- function(U, alpha = NULL, Cstar = NULL,
                               band = TRUE, conf = 0.95, B = 10000,
                               seed = 1, ...) {
  n <- length(U)
  if (n == 0) {
    message("No candidate T1Es detected; cannot produce QQ-plot.")
    return(invisible(NULL))
  }

  theoretical <- stats::ppoints(n)
  observed <- sort(U)

  # OLS best-fit line (computed first so the confidence band can follow it)
  fit_line <- stats::lm(observed ~ theoretical)
  a <- stats::coef(fit_line)[1]
  b <- stats::coef(fit_line)[2]

  # Simultaneous confidence band (follows best-fit line)
  # MC-calibrate pointwise level p* for Beta(i, n+1-i) bands such that
  # simultaneous coverage = conf, then transform bands through OLS fit.
  lower <- upper <- outside <- NULL
  if (band && n >= 2) {
    set.seed(seed)
    U_sim <- matrix(stats::runif(n * B), nrow = n, ncol = B)
    U_sim <- apply(U_sim, 2, sort)
    i <- seq_len(n)

    coverage_hat <- function(p) {
      tail <- (1 - p) / 2
      L <- stats::qbeta(tail, i, n + 1 - i)
      U <- stats::qbeta(1 - tail, i, n + 1 - i)
      mean(colSums(U_sim >= L & U_sim <= U) == n)
    }

    f <- function(p) coverage_hat(p) - conf
    if (f(conf) >= 0) {
      p_star <- conf
    } else {
      p_star <- stats::uniroot(f, lower = conf, upper = 0.9999, tol = 1e-4)$root
    }

    # Raw Unif(0,1) bands from Beta order statistics
    tail_star <- (1 - p_star) / 2
    lower_raw <- stats::qbeta(tail_star, i, n + 1 - i)
    upper_raw <- stats::qbeta(1 - tail_star, i, n + 1 - i)

    # Transform bands to follow the best-fit line
    lower <- a + b * lower_raw
    upper <- a + b * upper_raw

    outside_ord <- which(observed < lower | observed > upper)
    outside <- order(U)[outside_ord]
  }

  ## --- ggplot2 version ---
  main_title <- "Diagnostic QQ-Plot"
  if (!is.null(alpha) && !is.null(Cstar)) {
    main_title <- paste0("Diagnostic QQ-Plot at \u03b1 = ", alpha,
                         ", C* = ", round(Cstar, 4))
  } else if (!is.null(alpha)) {
    main_title <- paste0("Diagnostic QQ-Plot at \u03b1 = ", alpha)
  }

  df_plot <- data.frame(theoretical = theoretical, observed = observed)
  is_outside <- rep(FALSE, n)
  if (band && n >= 2) is_outside[outside_ord] <- TRUE
  df_plot$outside <- is_outside

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_point(ggplot2::aes(colour = outside), size = 1.5, shape = 16) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                                  guide = "none") +
    ggplot2::geom_abline(intercept = a, slope = b, colour = "red", linewidth = 1) +
    ggplot2::labs(x = "Theoretical Unif(0,1) Quantiles",
                  y = "Observed U_i Quantiles",
                  title = main_title) +
    ggplot2::theme_minimal()

  if (band && n >= 2) {
    df_band <- data.frame(theoretical = theoretical, lower = lower, upper = upper)
    p <- p +
      ggplot2::geom_line(data = df_band,
                          ggplot2::aes(x = theoretical, y = lower),
                          colour = "grey50", linetype = "dashed") +
      ggplot2::geom_line(data = df_band,
                          ggplot2::aes(x = theoretical, y = upper),
                          colour = "grey50", linetype = "dashed")
  }

  print(p)

  invisible(list(
    theoretical = theoretical,
    observed = observed,
    intercept = a,
    slope = b,
    lower = lower,
    upper = upper,
    outside = outside,
    plot = p
  ))
}


#' Adaptive calibration plot using C*(alpha)
#'
#' For each alpha on a grid in \code{[0, up]}, finds C*(alpha) such that
#' the observed proportion of contradictions (p <= alpha AND eJAB01 > C)
#' is closest to the target rate alpha/up.  Produces three plots:
#' \enumerate{
#'   \item Calibration curve: observed proportion vs alpha using the
#'         adaptive C*(alpha).  Should track the reference line (slope 1/up)
#'         closely by construction.
#'   \item C*(alpha) vs alpha: shows how the threshold changes with alpha
#'         (expected: decreasing).
#'   \item Diagnostic QQ-plot at the specified \code{alpha}, using
#'         C*(\code{alpha}).
#' }
#'
#' @param p Numeric vector of p-values (should span the full range used in
#'   estimation, not just significant results)
#' @param ejab Numeric vector of eJAB01 values (same length as \code{p})
#' @param up Upper bound for the alpha grid (default 0.1)
#' @param alpha Significance level for the diagnostic QQ-plot (default 0.05)
#' @param grid_range Length-2 numeric vector \code{c(lower, upper)} for the
#'   C search grid (default \code{c(0, 1)})
#' @param grid_n Number of C grid points (default 200)
#' @param n_alpha Number of alpha grid points (default 200)
#' @param n Numeric vector of sample sizes (same length as \code{p}); needed
#'   for the QQ-plot diagnostic.  If \code{NULL}, the QQ-plot is skipped.
#' @param q Numeric vector of parameter dimensions (same length as \code{p});
#'   needed for the QQ-plot diagnostic.  If \code{NULL}, the QQ-plot is skipped.
#' @param ... Additional arguments (currently unused, kept for compatibility)
#' @return Invisibly, a list with components:
#'   \describe{
#'     \item{alpha_grid}{Numeric vector of alpha values used.}
#'     \item{Cstar_alpha}{Numeric vector of C*(alpha) values.}
#'     \item{proportions}{Numeric vector of observed proportions at each alpha.}
#'     \item{Cstar_at_alpha}{C* at the specified \code{alpha}.}
#'   }
#' @export
calibration_plot <- function(p, ejab, up = 0.1, alpha = 0.05,
                              grid_range = c(0, 1), grid_n = 200,
                              n_alpha = 200, n = NULL, q = NULL, ...) {
  # --- Build grids ---
  alpha_grid <- seq(0, up, length.out = n_alpha + 1)[-1]  # skip alpha=0
  C_grid <- seq(grid_range[1], grid_range[2], length.out = grid_n)
  N <- sum(p < up)

  Cstar_alpha  <- numeric(length(alpha_grid))
  proportions  <- numeric(length(alpha_grid))

  # --- For each alpha, find C*(alpha) ---
  for (i in seq_along(alpha_grid)) {
    a <- alpha_grid[i]
    target <- a / up

    props <- vapply(C_grid, function(C) {
      sum(p <= a & ejab > C) / N
    }, numeric(1))

    best <- which.min((props - target)^2)
    Cstar_alpha[i] <- C_grid[best]
    proportions[i] <- props[best]
  }

  # C* at the specified alpha
  nearest <- which.min(abs(alpha_grid - alpha))
  Cstar_at_alpha <- Cstar_alpha[nearest]

  ## --- ggplot2 versions ---
  df_cal <- data.frame(alpha = alpha_grid, proportion = proportions)

  # Plot 1: Calibration curve
  p1 <- ggplot2::ggplot(df_cal, ggplot2::aes(x = alpha, y = proportion)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1 / up,
                          colour = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::xlim(0, up) +
    ggplot2::ylim(0, max(proportions, 1)) +
    ggplot2::labs(x = expression(alpha),
                  y = "Observed Proportion",
                  title = paste0("Calibration using C*(\u03b1)")) +
    ggplot2::theme_minimal()

  # Plot 2: C*(alpha) vs alpha
  df_cstar <- data.frame(alpha = alpha_grid, Cstar = Cstar_alpha)

  p2 <- ggplot2::ggplot(df_cstar, ggplot2::aes(x = alpha, y = Cstar)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey50", linetype = "dotted") +
    ggplot2::coord_cartesian(ylim = c(min(Cstar_alpha) * 0.9, max(Cstar_alpha) * 1.1)) +
    ggplot2::labs(x = expression(alpha),
                  y = expression(C^"*" * (alpha)),
                  title = paste0("C*(\u03b1) vs \u03b1")) +
    ggplot2::theme_minimal()

  # Plot 3: Diagnostic QQ-plot at specified alpha
  p3 <- NULL
  if (!is.null(n) && !is.null(q)) {
    idx <- which(p < alpha & ejab > Cstar_at_alpha)
    if (length(idx) > 0) {
      U <- diagnostic_U(p[idx], n[idx], q[idx], alpha, Cstar_at_alpha)
      qq_result <- diagnostic_qqplot(U, alpha = alpha, Cstar = Cstar_at_alpha)
      p3 <- qq_result$plot
    } else {
      message("No candidates at alpha = ", alpha, " with C* = ",
              round(Cstar_at_alpha, 4))
    }
  }

  print(p1)
  print(p2)

  invisible(list(alpha_grid = alpha_grid,
                 Cstar_alpha = Cstar_alpha,
                 proportions = proportions,
                 Cstar_at_alpha = Cstar_at_alpha,
                 plots = list(calibration = p1, cstar = p2, qqplot = p3)))
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

  # Calibration plot (includes QQ-plot as 3rd panel)
  cal <- NULL
  if (plot) {
    cal <- calibration_plot(df$p, ejab_vals, up = up, alpha = alpha,
                            n = df$n, q = df$q)
  }

  # Compute diagnostics for candidates
  Cstar_alpha <- if (!is.null(cal)) cal$Cstar_at_alpha else Cfit$Cstar
  U <- NULL
  if (length(idx) > 0) {
    U <- diagnostic_U(df$p[idx], df$n[idx], df$q[idx], alpha, Cfit$Cstar)
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
