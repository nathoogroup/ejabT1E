\name{ejab01}
\alias{ejab01}
\title{Compute eJAB01}
\description{
Computes the approximate objective Bayes factor eJAB01 for each NHST result.
}
\usage{
ejab01(p, n, q)
}
\arguments{
\item{p}{Numeric vector of p-values (0 < p < 1).}
\item{n}{Numeric vector of sample sizes (per test, must be > 1).}
\item{q}{Numeric vector of parameter dimensions (per test, >= 1).}
}
\value{
Numeric vector of eJAB01 values.
}
\examples{
ejab01(0.03, 50, 1)
ejab01(c(0.01, 0.04), c(100, 200), c(1, 2))
}
