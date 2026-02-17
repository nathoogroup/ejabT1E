\name{diagnostic_qqplot}
\alias{diagnostic_qqplot}
\title{Produce Diagnostic QQ-plot}
\description{
Produces a QQ-plot of the diagnostic U_i values against
Unif(0,1) theoretical quantiles. Linearity indicates that
the left-tail uniformity assumption holds. A best-fit line
(via OLS) is overlaid instead of the 45-degree line, since
C* is estimated. Optionally adds a simultaneous confidence
band based on Beta order statistics, transformed to follow the
best-fit line rather than the 45-degree line.
}
\usage{
diagnostic_qqplot(U, alpha = NULL, Cstar = NULL,
                  band = TRUE, conf = 0.95, B = 10000,
                  seed = 1, ...)
}
\arguments{
\item{U}{Numeric vector of diagnostic U_i values.}
\item{alpha}{Significance level used for detection (displayed in title).
  If NULL, omitted from title.}
\item{Cstar}{C* value used for detection (displayed in title).
  If NULL, omitted from title.}
\item{band}{Logical; add simultaneous confidence band? (default TRUE).}
\item{conf}{Confidence level for band (default 0.95).}
\item{B}{Number of Monte Carlo simulations for band calibration (default 10000).}
\item{seed}{Random seed for reproducibility (default 1).}
\item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
Invisibly, a list with components:
\describe{
  \item{theoretical}{Theoretical quantiles.}
  \item{observed}{Observed (sorted) U values.}
  \item{intercept}{OLS intercept.}
  \item{slope}{OLS slope.}
  \item{lower}{Lower band (if band=TRUE).}
  \item{upper}{Upper band (if band=TRUE).}
  \item{outside}{Indices of observations outside band (if band=TRUE).}
}
}
