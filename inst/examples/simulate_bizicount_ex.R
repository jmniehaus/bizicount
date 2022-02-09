## SETUP
set.seed(123)
n = 150

# define a function to simulate from a gaussian copula
# first margin is zero-inflated negative binomial (zinb)
# second margin is zero-inflated poisson (zip)
# Note: marginal distributions are hard-coded in function, including
# inverse dispersion parameter for zinb.
gen = function(n,
               b1,
               b2,
               g1,
               g2,
               dep) {

     k1 = length(b1)
     k2 = length(b2)

     X1 = cbind(1, matrix(rbinom(n * (k1 - 1), 1, .5), ncol = k1 - 1))
     X2 = cbind(1, matrix(rexp(n * (k2 - 1), 3), ncol = k2 - 1))

     lam1 = exp(X1 %*% b1)
     lam2 = exp(X2 %*% b2)

     Z1 = cbind(1, matrix(runif(n * (k1 - 1), -1, 1), ncol = k1 - 1))
     Z2 = cbind(1, matrix(rnorm(n * (k2 - 1)), ncol = k2 - 1))

     psi1 = plogis(Z1 %*% g1)
     psi2 = plogis(Z2 %*% g2)

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

# simulate from fitted model
sims = simulate(mod, nsim = 150)


# input sims to DHARMa for diagnostics
# margin 1
d1 = DHARMa::createDHARMa(
     simulatedResponse = sims[[1]],
     observedResponse = dat$y1,
     fittedPredictedResponse = fitted(mod)[,1],
     integerResponse = TRUE,
     method = "PIT"
)

# margin 2
d2 = DHARMa::createDHARMa(
     simulatedResponse = sims[[2]],
     observedResponse = dat$y2,
     fittedPredictedResponse = fitted(mod)[,2],
     integerResponse = TRUE,
     method = "PIT"
)

# test each margin
DHARMa::testResiduals(d1)
DHARMa::testResiduals(d2)


