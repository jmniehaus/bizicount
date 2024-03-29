% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zi_dist_fns.R
\name{zip}
\alias{zip}
\alias{dzip}
\alias{rzip}
\alias{pzip}
\alias{qzip}
\title{The zero-inflated Poisson (ZIP) distribution}
\usage{
dzip(x, lambda, psi, log = FALSE, recycle = FALSE)

rzip(n, lambda, psi, recycle = FALSE)

pzip(q, lambda, psi, lower.tail = TRUE, log.p = FALSE, recycle = FALSE)

qzip(p, lambda, psi, lower.tail = TRUE, log.p = FALSE, recycle = FALSE)
}
\arguments{
\item{x, q}{Vector of quantiles at which to evaluate the PMF and CDF,
respectively. Should be non-negative integers.}

\item{lambda}{Vector of means for the count portion of the zero-inflated
Poisson distribution. Should be non-negative. NOTE: This is \emph{not} the mean
of the zero-inflated Poisson distribution; it is the mean of the Poisson
component of the mixture distribution. See 'Details.'}

\item{psi}{Vector of zero-inflation probabilities.}

\item{log, log.p}{Logical indicating whether probabilities should be returned
on log scale (for \code{dzip} and \code{pzip}), or are supplied on log-scale (for \code{qzip}).}

\item{recycle}{Logical indicating whether to permit arbitrary recycling of
arguments with unequal length. See 'Details' and 'Examples.'}

\item{n}{Number of realizations from the distribution to generate}

\item{lower.tail}{Logical indicating whether probabilities should be
\eqn{Pr(X \le x)} or \eqn{Pr(X > x)}}

\item{p}{Vector of probabilities at which to evaluate the quantile function.}
}
\value{
\code{dzip} returns the mass function evaluated at \code{x},
\code{pzip} returns the CDF evaluated at \code{q}, \code{qzip} returns the quantile
function evaluated at \code{p}, and \code{rzip} returns random variates with the
specified parameters.
}
\description{
These functions are used to evaluate the zero-inflated Poisson
distribution's probability mass function (PMF), cumulative distribution
function (CDF), and quantile function (inverse CDF), as well as generate
random realizations from the ZIP distribution.
}
\details{
The zero inflated Poisson distribution is a mixture of a Poisson
and a degenerate point-mass at 0. It has the form
\deqn{\psi + (1-\psi)(\lambda^x e^-\lambda)/x!}, with mean
\eqn{(1-\psi)\lambda}. Thus, the parameter \code{lambda} above is the mean of
the Poisson distribution that forms part of the zero-inflated
distribution, \emph{not} the mean of the ZIP distribution.

\code{recycle} -- If \code{FALSE} (default), all arguments must have identical
length, there can be two unique lengths for the arguments, provided that
one of those lengths is 1. For example, \code{lambda = c(1,2,3)} and \code{psi=.5}
is acceptable because there are two unique lengths, and one of them is
length 1. However, \code{lambda=c(1,2,3)} and \code{psi=c(.5,.2)} would fail, as
there are two distinct lengths, none of which is 1. If \verb{TRUE,} no
additional checks (beyond those in base R's functions) are made to ensure
that the argument vectors have the same length.
}
\examples{
# Unequal lengths, but one of them is length 1, others are same length (3).
# No error.

x = c(1,2,3)
lambda = c(3,4,5)
psi = .1

dzip(x, lambda, psi)


# unequal lengths, at least one of them is not length 1,
# error

\dontrun{

x = c(1,2,3)
lambda = c(3,4)
psi = .1

dzip(x, lambda, psi)

}

# unequal lengths, at least one of them is not length 1.
# but set recycle = T to permit arbitrary recycling.

x = c(1,2,3)
lambda = c(3,4)
psi = .1

dzip(x, lambda, psi, recycle=TRUE)
}
\references{
Lambert, Diane. "Zero-inflated Poisson regression, with an
application to defects in manufacturing." Technometrics 34.1 (1992): 1-14.
}
\author{
John Niehaus
}
