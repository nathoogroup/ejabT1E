\name{objective_C}
\alias{objective_C}
\title{Closed-form Objective Function for C}
\description{
Evaluates the integrated squared difference between the empirical
contradiction rate and the theoretical Type I error rate alpha/u_p,
using the closed-form expression derived from the step-function
structure of P_hat.
}
\usage{
objective_C(C, p, ejab, up)
}
\arguments{
\item{C}{Candidate threshold value.}
\item{p}{Numeric vector of all p-values.}
\item{ejab}{Numeric vector of all eJAB01 values.}
\item{up}{Upper bound for p-value filtering.}
}
\details{
The objective is:
\deqn{\sum_{j=1}^{J} (p_{(j+1)} - p_{(j)}) \hat{P}(p_{(j)}, C)^2
  - \frac{1}{u_p N} \sum_{j=1}^{J} \xi_j^C (u_p^2 - p_{(j)}^2) + \frac{u_p}{3}}
where \eqn{p_{(J+1)} = u_p}, \eqn{\xi_j^C} counts results at p-value \eqn{p_{(j)}}
with eJAB > C, and \eqn{\hat{P}(p_{(j)}, C) = \frac{1}{N}\sum_{k=1}^{j}\xi_k^C}.
}
\value{
Scalar objective function value.
}
\examples{
set.seed(1)
p <- runif(100, 0, 0.05)
e <- ejab01(p, rep(50, 100), rep(1, 100))
objective_C(1.0, p, e, up = 0.05)
}
