set.seed(123)
n = 500
u = rpois(n, 3)
y1 = rzip(n, 12, .2) + u
y2 = rpois(n, 8) + u

# Single parameter test, covariates can be added though.
uni1 = glm(y1 ~ 1, family = poisson())
uni2 = glm(y2 ~ 1, family = poisson())

biv = bizicount(y1~1, y2~1, margins = c("pois", "pois"), keep = TRUE)

zi_test(uni1)
zi_test(uni2)

zi_test(biv)
