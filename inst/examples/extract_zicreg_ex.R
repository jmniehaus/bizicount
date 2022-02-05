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

mod = zic.reg(y ~ x | z, data = dat)


### Output to table with texreg

# extract information

tr_obj_se = texreg::extract(mod)
tr_obj_ci = texreg::extract(mod, CI = .95)

# output to latex, single table

texreg::texreg(list(tr_obj_se, tr_obj_ci))

# output to plain text, multiple tables

texreg::screenreg(tr_obj_se)
texreg::screenreg(tr_obj_ci)
