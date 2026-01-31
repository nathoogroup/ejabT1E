\name{posterior_t1e}
\alias{posterior_t1e}
\title{Compute Posterior Probability of Type I Error}
\description{
Treats C* as the posterior mode under a uniform prior on the grid.
Converts objective function values to a posterior distribution over C,
then for each result computes P(T1E) by summing posterior mass over
all C values satisfying the contradiction criterion.
}
\usage{
posterior_t1e(p, ejab, alpha, all_objectives, grid)
}
\arguments{
\item{p}{Numeric vector of p-values (for results of interest).}
\item{ejab}{Numeric vector of eJAB01 values.}
\item{alpha}{Significance level.}
\item{all_objectives}{Numeric vector of objective function values over the grid.}
\item{grid}{Numeric vector of C values corresponding to all_objectives.}
}
\details{
The posterior over C is obtained by:
\deqn{\pi(C_k | \text{data}) \propto \exp\{-\text{objective}(C_k)\}}
assuming a uniform prior over the grid. For each result i, the
posterior probability of being a Type I error is:
\deqn{P(\text{T1E}_i) = \sum_{k: p_i \le \alpha,\; \text{eJAB}_i > C_k} \pi(C_k | \text{data})}
}
\value{
Numeric vector of posterior T1E probabilities, one per result.
}
