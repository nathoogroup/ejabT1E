\name{ejab_pipeline}
\alias{ejab_pipeline}
\title{Full eJAB Type I Error Detection Pipeline}
\description{
Runs the complete analysis: computes eJAB01 values, estimates C*,
detects candidate Type I errors, and produces diagnostics.
}
\usage{
ejab_pipeline(df, up = 0.05, alpha = 0.05,
              grid_range = c(1/3, 3), grid_n = 200, plot = TRUE)
}
\arguments{
\item{df}{Data frame with columns: \code{p} (p-values), \code{n}
  (sample sizes, > 1), \code{q} (test dimensions). An optional \code{ID}
  column is preserved in output.}
\item{up}{Upper p-value bound (default 0.05). Only results with
  p < up are used for estimating C*.}
\item{alpha}{Significance level for T1E detection (default 0.05).
  Must satisfy alpha <= up.}
\item{grid_range}{Length-2 numeric vector \code{c(lower, upper)} for the
  C* grid (default \code{c(1/3, 3)}). Use \code{c(1, 1)} to fix C* = 1.}
\item{grid_n}{Number of grid points (default 200).}
\item{plot}{Logical; produce calibration plot and diagnostic QQ-plot?
  (default TRUE).}
}
\value{
A list with components:
\describe{
  \item{Cstar}{Estimated threshold C*.}
  \item{objective}{Minimized objective function value.}
  \item{candidates}{Data frame of candidate T1Es with columns from input
    plus \code{ejab}. NULL if none detected.}
  \item{U}{Numeric vector of diagnostic U_i values for candidates.}
  \item{ejab}{Numeric vector of eJAB01 values for all input results.}
}
}
\examples{
set.seed(42)
df <- data.frame(
  ID = 1:100,
  p = runif(100, 0, 0.05),
  n = sample(20:200, 100, replace = TRUE),
  q = rep(1, 100)
)
result <- ejab_pipeline(df, up = 0.05, alpha = 0.05, plot = FALSE)
result$Cstar
head(result$candidates)
}
