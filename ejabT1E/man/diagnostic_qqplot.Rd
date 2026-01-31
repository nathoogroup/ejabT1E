\name{diagnostic_qqplot}
\alias{diagnostic_qqplot}
\title{Produce Diagnostic QQ-plot}
\description{
Produces a QQ-plot of the diagnostic U_i values against
Unif(0,1) theoretical quantiles. Linearity indicates that
the left-tail uniformity assumption holds.
}
\usage{
diagnostic_qqplot(U, ...)
}
\arguments{
\item{U}{Numeric vector of diagnostic U_i values.}
\item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}
\value{
Invisibly, a list with components \code{theoretical} and \code{observed}.
}
