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


#' Estimate C* for different alpha values
#'
#' For each alpha value on a grid, finds the C value that minimizes the squared
#' difference between the observed proportion of contradictions and the target
#' proportion alpha/up.
#'
#' @param p Numeric vector of p-values
#' @param ejab Numeric vector of eJAB01 values (same length as \code{p})
#' @param up Upper bound for the alpha grid (default 0.1)
#' @param grid_range Length-2 numeric vector \code{c(lower, upper)} for the
#'   C grid (default \code{c(0, 1)})
#' @param grid_n Number of C grid points (default 200)
#' @param n_alpha Number of alpha grid points (default 200)
#' @return A list with components:
#'   \describe{
#'     \item{alpha_grid}{Numeric vector of alpha values}
#'     \item{Cstar_alpha}{Numeric vector of C* values, one per alpha}
#'     \item{proportions}{Numeric vector of observed proportions at each alpha}
#'   }
#' @export
estimate_Cstar_alpha <- function(p, ejab, up = 0.1,
                                  grid_range = c(0, 1), grid_n = 200,
                                  n_alpha = 200) {
  alpha_grid <- seq(0, up, length.out = n_alpha + 1)[-1]
  C_grid     <- seq(grid_range[1], grid_range[2], length.out = grid_n)
  N          <- sum(p < up)

  Cstar_alpha <- numeric(length(alpha_grid))
  proportions <- numeric(length(alpha_grid))

  for (i in seq_along(alpha_grid)) {
    a      <- alpha_grid[i]
    target <- a / up
    props  <- vapply(C_grid, function(C) sum(p <= a & ejab > C) / N, numeric(1))
    best   <- which.min((props - target)^2)
    Cstar_alpha[i] <- C_grid[best]
    proportions[i] <- props[best]
  }

  list(alpha_grid = alpha_grid,
       Cstar_alpha = Cstar_alpha,
       proportions = proportions)
}
