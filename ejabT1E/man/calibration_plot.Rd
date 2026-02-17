\name{calibration_plot}
\alias{calibration_plot}
\title{Adaptive calibration plot using C*(alpha)}
\description{
For each alpha on a grid in [0, up], finds C*(alpha) such that the
observed proportion of contradictions (p <= alpha AND eJAB01 > C) is
closest to the target rate alpha/up.  Produces three plots:
(1) calibration curve using adaptive C*(alpha), which should track
the reference line (slope 1/up) closely;
(2) C*(alpha) vs alpha, showing how the threshold changes (expected:
decreasing);
(3) diagnostic QQ-plot at the specified \code{alpha} using C*(alpha).
}
\usage{
calibration_plot(p, ejab, up = 0.1, alpha = 0.05,
                 grid_range = c(0, 1), grid_n = 200,
                 n_alpha = 200, n = NULL, q = NULL, ...)
}
\arguments{
\item{p}{Numeric vector of p-values (should span the full range used in
  estimation, not just significant results).}
\item{ejab}{Numeric vector of eJAB01 values (same length as \code{p}).}
\item{up}{Upper bound for the alpha grid (default 0.1).}
\item{alpha}{Significance level for the diagnostic QQ-plot (default 0.05).}
\item{grid_range}{Length-2 numeric vector \code{c(lower, upper)} for the
  C search grid (default \code{c(0, 1)}).}
\item{grid_n}{Number of C grid points (default 200).}
\item{n_alpha}{Number of alpha grid points (default 200).}
\item{n}{Numeric vector of sample sizes (same length as \code{p}); needed
  for the QQ-plot diagnostic.  If \code{NULL}, the QQ-plot is skipped.}
\item{q}{Numeric vector of parameter dimensions (same length as \code{p});
  needed for the QQ-plot diagnostic.  If \code{NULL}, the QQ-plot is skipped.}
\item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
Invisibly, a list with components:
\describe{
  \item{alpha_grid}{Numeric vector of alpha values used.}
  \item{Cstar_alpha}{Numeric vector of C*(alpha) values.}
  \item{proportions}{Numeric vector of observed proportions at each alpha.}
  \item{Cstar_at_alpha}{C* at the specified \code{alpha}.}
}
}
