# Simulate some zip data
n=1000
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


### Predict on observed/training data
# predict conditional mean (fitted values)
predict(mod, type = "mean")

# predict probabilty Y = y
probs_pred_obs = predict(mod, type = "prob")

# predict mean of count distribution (lambda)
lambda_pred_obs = predict(mod, type = "lambda")

# mse predicted vs true lambda values
mean((lam - lambda_pred_obs)**2)

# predict zero inflation probability (psi)
psi_pred_obs = predict(mod, type = "psi")

# MSE predicted vs true zero-inflation probabilities
mean((psi-psi_pred_obs)**2)


### Predict on test data
# simulate some test data

x = cbind(1, rnorm(n, mean = -0.5, sd = 1.25))
z = cbind(1, rbeta(n, 6, 12))
y = rzip(n, lambda = exp(x %*% coef(mod)[1:2]), psi = plogis(z %*% coef(mod)[3:4]))
dat_new = cbind.data.frame(x = x[,-1], z = z[,-1], y = y)

# predict conditional mean
mean_new = predict(mod, type = "mean", newdata = dat_new)
mean((y - mean_new)**2)

# predict probability of Y = y
probs_new = predict(mod, type = "prob", newdata = dat_new, y.new = y)

# predict lambda
lambda_new = predict(mod, type = "lambda", newdata = dat_new)

# predict zero inflation probability
psi_new = predict(mod, type = "psi", newdata = dat_new)




