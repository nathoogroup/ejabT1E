\name{calibration_plot}
\alias{calibration_plot}
\title{Calibration plot (Figure 3a style)}
\description{
Multi-panel calibration plot showing observed proportion of contradictions
(p <= alpha AND eJAB01 > C) against alpha for a grid of C values. If the
method is well-calibrated, the curve should follow the reference line
with slope 1/up. Use \code{grid_range = c(C, C)} for a single-panel plot
at a specific C value.
}
\usage{
calibration_plot(p, ejab, up = 0.1, grid_range = c(0, 1),
                 grid_n = 18, n_alpha = 200, ...)
}
\arguments{
\item{p}{Numeric vector of p-values (should span the full range used in
  estimation, not just significant results).}
\item{ejab}{Numeric vector of eJAB01 values (same length as \code{p}).}
\item{up}{Upper bound for the alpha grid (default 0.1).}
\item{grid_range}{Length-2 numeric vector \code{c(lower, upper)} for the
  C grid (default \code{c(0, 1)}). Use \code{c(C, C)} for a single panel.}
\item{grid_n}{Number of C grid points (default 18). Panels are arranged
  in a 3x3 grid; more than 9 panels will span multiple pages.}
\item{n_alpha}{Number of alpha grid points per panel (default 200).}
\item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
Invisibly, a list with components \code{C_grid} (numeric vector of
C values), \code{alpha_grid} (numeric vector of alpha values), and
\code{proportions} (matrix with one row per C value and one column per
alpha value).
}
