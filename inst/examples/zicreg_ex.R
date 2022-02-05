## ZIP example
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

# estimate zip model using NLM, no data.frame

mod = zic.reg(y ~ x[,-1] | z[,-1])

# same model, with dataframe

mod = zic.reg(y ~ x | z, data = dat)


# estimate zip using NLM, adjust stepmax via ... param

mod = zic.reg(y ~ x[,-1] | z[,-1], stepmax = .5)


# estimate zip using optim

mod = zic.reg(y ~ x[,-1] | z[,-1], optimizer = "optim")


# pass different method, reltol to optim using ... param

mod = zic.reg(y ~ x[,-1] | z[,-1],
        optimizer = "optim",
        method = "Nelder-Mead",
        control = list(reltol = 1e-10)
        )

# No formula, specify design matrices and offsets.
zic.reg(y=y, X=x, z=z)



## ZINB example
# simulate zinb data

disp = .5
y = bizicount::rzinb(n, psi = psi, size = disp, mu=lam)


# zinb model, use keep = TRUE for post-estimation methods

mod = zic.reg(y ~ x[,-1] | z[,-1], dist = "n", keep = TRUE)

print(mod)
summary(mod)

