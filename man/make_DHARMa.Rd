% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bizicount-methods.R
\name{make_DHARMa}
\alias{make_DHARMa}
\title{DHARMa-class objects from bizicount models}
\usage{
make_DHARMa(object, nsim = 250, seed = 123, method = "PIT")
}
\arguments{
\item{object}{A \code{\link{bizicount-class}} object, as returned by \link{bizicount}.}

\item{nsim}{Number of simulated responses from the fitted model to use for diagnostics.}

\item{seed}{Random seed for simulating from fitted model.}

\item{method}{See \code{\link[DHARMa]{createDHARMa}}.}
}
\value{
A list of \code{\link[DHARMa]{DHARMa}} objects.
}
\description{
A wrapper around the \code{\link[DHARMa]{DHARMa}} package's
\code{\link[DHARMa]{createDHARMa}} function. Creates a list of DHARMa
objects, one for each margin of a \code{\link{bizicount-class}} object, using
simulated responses from \code{\link{simulate.bizicount}}.
}
\note{
This is merely a wrapper around the \code{\link[DHARMa]{createDHARMa}}
function to conveniently get DHARMa objects for each margin of a bizicount
model.
}
\examples{
## SETUP
set.seed(123)
n = 100

# define a function to simulate from a gaussian copula
# first margin is zero-inflated negative binomial (zinb)
# second margin is zero-inflated poisson (zip)
# Note: marginal distributions are hard-coded in function, including
# inverse dispersion parameter for zinb.
gen = function(n, b1, b2, g1, g2, dep) {

     k1 = length(b1)
     k2 = length(b2)

     X1 = cbind(1, matrix(rbinom(n * (k1 - 1), 1, .5), ncol = k1 - 1))
     X2 = cbind(1, matrix(rexp(n * (k2 - 1), 3), ncol = k2 - 1))

     lam1 = exp(X1 \%*\% b1)
     lam2 = exp(X2 \%*\% b2)

     Z1 = cbind(1, matrix(runif(n * (k1 - 1), -1, 1), ncol = k1 - 1))
     Z2 = cbind(1, matrix(rnorm(n * (k2 - 1)), ncol = k2 - 1))

     psi1 = plogis(Z1 \%*\% g1)
     psi2 = plogis(Z2 \%*\% g2)

     norm_vars = MASS::mvrnorm(
          n,
          mu = c(0, 0),
          Sigma = matrix(c(1, dep, dep, 1), ncol =2)
     )

     U = pnorm(norm_vars)

     y1 =  qzinb(U[, 1],
                 mu = lam1,
                 psi = psi1,
                 size = .3)
     y2 =  qzip(U[, 2],
                lambda = lam2,
                psi = psi2)

     dat = data.frame(
          X1 = X1[, -1],
          X2 = X2[, -1],
          Z1 = Z1[, -1],
          Z2 = Z2[, -1],
          y1,
          y2,
          lam1,
          lam2,
          psi1,
          psi2
     )
     return(dat)
}


# define parameters
b1 = c(1, -2, 3)
b2 = c(-1, 3, 1)
g1 = c(2, -1.5, 2)
g2 = c(-1, -3.75, 1.25)
rho = .5


# generate data
dat = gen(n, b1, b2, g1, g2, rho)
f1 = y1 ~ X1.1 + X1.2 | Z1.1 + Z1.2
f2 = y2 ~ X2.1 + X2.2 | Z2.1 + Z2.2

## END SETUP




# estimate model

mod = bizicount(f1, f2, dat, cop = "g", margins = c("zinb", "zip"), keep=TRUE)


# diagnose model with DHARMa
# see end for simulate.bizicount example.

dharm = make_DHARMa(mod, nsim = 100)

lapply(dharm, DHARMa::testResiduals)
}
\references{
Florian Hartig (2022). DHARMa: Residual Diagnostics for
Hierarchical (Multi-Level / Mixed) Regression Models. R package version
0.4.5. https://CRAN.R-project.org/package=DHARMa
}
\seealso{
\code{\link[DHARMa]{DHARMa}}, \code{\link[DHARMa]{createDHARMa}},
\code{\link{simulate.bizicount}}
}
\author{
John Niehaus
}
