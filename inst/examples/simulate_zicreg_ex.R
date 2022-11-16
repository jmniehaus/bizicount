# Simulate some zip data
n=300
x = cbind(1, rnorm(n))
z = cbind(1, rbeta(n, 4, 8))
b = c(1, 2.2)
g = c(-1, 1.7)
lam = exp(x %*% b)
psi = plogis(z %*% g)

y = bizicount::rzip(n, lambda = lam, psi=psi)
dat = cbind.data.frame(x = x[,-1], z = z[,-1], y = y)

# estimate model

mod = zic.reg(y ~ x | z, data = dat, keep = TRUE)

# simulate from fit for use in dharma
sims = simulate(mod)

### Make dharma object

dharm = DHARMa::createDHARMa(
     simulatedResponse = sims,
     observedResponse = y,
     fittedPredictedResponse = fitted(mod),
     integerResponse = TRUE,
     method = "PIT"
)

### Plot the DHARMa object, do other diagnostics
plot(dharm)
DHARMa::testResiduals(dharm)
