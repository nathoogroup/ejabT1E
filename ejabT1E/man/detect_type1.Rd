\name{detect_type1}
\alias{detect_type1}
\title{Detect Candidate Type I Errors}
\description{
Identifies results that are Bayes/NHST contradictions: p-value <= alpha
(NHST rejects) and eJAB01 > Cstar (Bayes factor supports H0).
}
\usage{
detect_type1(p, ejab, alpha, Cstar)
}
\arguments{
\item{p}{Numeric vector of p-values.}
\item{ejab}{Numeric vector of eJAB01 values.}
\item{alpha}{Significance level (must be <= up used in estimation).}
\item{Cstar}{Estimated threshold.}
}
\value{
Integer vector of indices of candidate Type I errors.
}
\examples{
set.seed(1)
p <- runif(100, 0, 0.05)
e <- ejab01(p, rep(50, 100), rep(1, 100))
idx <- detect_type1(p, e, alpha = 0.05, Cstar = 1.0)
length(idx)
}
