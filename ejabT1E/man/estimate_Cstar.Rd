\name{estimate_Cstar}
\alias{estimate_Cstar}
\title{Grid Search for C*}
\description{
Finds the value of C on a grid that minimizes the objective function.
The default grid is \code{seq(1/3, 3, length.out = 200)}.
}
\usage{
estimate_Cstar(p, ejab, up = 0.05, grid = seq(1/3, 3, length.out = 200))
}
\arguments{
\item{p}{Numeric vector of p-values.}
\item{ejab}{Numeric vector of eJAB01 values.}
\item{up}{Upper bound (default 0.05).}
\item{grid}{Numeric vector of candidate C values
  (default: \code{seq(1/3, 3, length.out = 200)}).}
}
\value{
A list with components:
\describe{
  \item{Cstar}{The C value minimizing the objective.}
  \item{objective}{The minimized objective function value.}
  \item{all_objectives}{Numeric vector of objective values over the full grid.}
  \item{grid}{The grid of C values used.}
}
}
\examples{
set.seed(1)
p <- runif(200, 0, 0.05)
e <- ejab01(p, rep(50, 200), rep(1, 200))
fit <- estimate_Cstar(p, e, up = 0.05)
fit$Cstar
}
