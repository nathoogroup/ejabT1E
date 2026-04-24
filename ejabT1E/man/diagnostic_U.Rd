\name{diagnostic_U}
\alias{diagnostic_U}
\title{Compute Diagnostic U_i Values}
\description{
For each candidate Type I error, computes a diagnostic U_i that should
follow Unif(0,1) if the left-tail uniformity assumption holds and the
selected results are true Type I errors.
}
\usage{
diagnostic_U(p, n, q, alpha, Cstar)
}
\arguments{
\item{p}{Numeric vector of p-values (for detected T1Es only).}
\item{n}{Numeric vector of sample sizes (per test, must be > 1).}
\item{q}{Numeric vector of dimensions (per test).}
\item{alpha}{Significance level.}
\item{Cstar}{Estimated threshold.}
}
\details{
The diagnostic is defined as:
\deqn{U_i = \frac{p_i - d_i}{\alpha - d_i}}
where
\deqn{d_i = 1 - F_{\chi^2_{q_i}}\left(\frac{2 n_i^{1/q_i}}{n_i^{1/q_i} - 1}
  \log\frac{\sqrt{n_i}}{C^*}\right)}
}
\value{
Numeric vector of U_i values.
}
