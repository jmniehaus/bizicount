pks = c("copula", "pbivnorm", "Formula",
        "numDeriv", "compiler", "microbenchmark", "texreg",
        "rlang")
needed = setdiff(pks, rownames(installed.packages()))
if(length(needed) > 0){
  options(Ncpus=(parallel::detectCores() -1))
  install.packages(needed)
  }
invisible(sapply(pks, library, character.only=T))

scripts = list.files(path="scripts", pattern="[.]R$", full.names=T)
scripts = grep("replication|plots", scripts, invert=T, value=T)
invisible(sapply(scripts, source))


### testing bivariate copula regression



#
set.seed(123)
n=500

#count coefs



u = rpois(n,2)#cbind(1, rnorm(n)) %*% c(1,-1.25)
y1 = rpois(n, 3) + u
y2 = rpois(n, 7) + u
X1 = cbind(1, rnorm(n))
X2 = cbind(1, runif(n, -1, 1))

b1 = c(1.44, 1.1)
lam1 = exp(X1 %*% b1)
lam1 = exp(1.44*X1[,2])

b2 = c(1.83, 2.63)
lam2 = exp(X2 %*% b2)
lam2 = exp(2.63*X2[,2])

#zi coefs
Z1 = cbind(1, rnorm(n))
Z2 = cbind(1, rt(n,5))

g1 = c(-2.64, .42)
psi1 = plogis(Z1 %*% g1)

g2 = c(-1.65, -2.75)
psi2 = plogis(Z2 %*% g2)

u = rzip(n, 2, .1)
y1 = rzip(n, lambda=5, .15) + u
y2 = rzip(n, lambda=3, .25) + u
y1 = rzip(n, lambda=lam1, psi1) +u
y2 = rzip(n, lambda=lam2, psi2) + u

y1 = rpois(n, 5)
y2 = rpois(n, 3)

y1 = rpois(n, lam1) + u
y2 = rpois(n, lam2) + u

gen = function(n,
               b1,
               b2,
               g1,
               g2,
               dep){

  k1 = length(b1)
  k2 = length(b2)
  X1 = cbind(1, matrix(rbinom(n*(k1-1), 1, .5), ncol = k1-1))
  X2 = cbind(1, matrix(rexp(n*(k2-1), 3), ncol = k2-1))
  lam1 = exp(X1%*%b1)
  lam2 = exp(X2%*%b2)
  Z1 = cbind(1, matrix(runif(n*(k1-1), -1, 1), ncol=k1-1))
  Z2 = cbind(1, matrix(rnorm(n*(k2-1)), ncol=k2-1))
  psi1 = plogis(Z1 %*% g1)
  psi2 = plogis(Z2 %*% g2)

  norm_vars = MASS::mvrnorm(n, mu = c(0,0), Sigma=matrix(c(1,dep,dep,1), ncol=2))
  U = pnorm(norm_vars)
  y1 =  qzinb(U[,1], mu=lam1, psi=psi1, size=.3)
  y2 =  qzinb(U[,2], mu=lam2, psi=psi2, size=.15)

  dat = data.frame(X1 = X1[,-1], X2 = X2[,-1], Z1=Z1[,-1], Z2=Z2[,-1], y1, y2, lam1, lam2, psi1, psi2)
  return(dat)
}
m = zic.reg(f1, data=data, dist="nb", keep=T)
d = createDHARMa(simulatedResponse= simulate(m),
             observedResponse = data$y1,
             fittedPredictedResponse = fitted(m),
             integerResponse = T,
             seed = 123,
             method = "PIT")
b1 = c(1, -2, 3)
b2 = c(-1, 3, 1)
g1 = c(-2, -1.5, 4)
g2 = c(-1, -3.75, 1.25)

data = gen(1500, b1, b2, g1, g2, .5)
f1 = y1 ~ X1 | Z1
f2 = y2 ~ X2 | Z2

y1 = rzinb(500, size=.3, psi=.15, mu=8)
y2 = rzinb(500, size=.15, psi=.3, mu=12)

mod = bizicount(y1~1|1,
                y2~1|1,
                margins=c("zinb", "zinb"))

mod = bizicount(f1, f2,
                margins = c("zip", "zinb"),
                cop = "g",
                data=data
                )

mod = bizicount(y1~X1[,2]|Z1[,2],
                y2~X2[,2]|Z2[,2],
                margins = c("zip", "zip")
                )

head(exp(1.44*mean(X1[,2]))*1.44)
head(exp(1.253 + .77*mean(X1[,2]))*.77)


ours = bizicount(y1~1|1,
                  y2~1|1,
                 margins = c("zip", "zip"),
          cop="frank",
          keep=T)

ours = bizicount(y1~1,
                 y2~1,
                 cop="frank")

gj = gjrm(list(y1~X1[,2], y2~X2[,2]), BivD = "F", Model="B", margins=c("PO", "PO"))
rbvpois
cbind(
coef(pscl::zeroinfl(formula = y1 ~ X1[, -1] | Z1[, -1])),
coef(zic.reg(y1~X1[,-1]|Z1[,-1])),
c(b1,g1)
)

cbind(
coef(bizicount(y[,-1]~X1[,-1]|Z1[,-1],
          y[,-2]~X2[,-1]|Z2[,-1],
          cop = "gaus",
          margins = c("zip", "zip"))),
c(b1,b2,g1,g2,NA)
)

dat = corrvar2(n=n,
               k_pois=2,
               lam = cbind(lam1, lam2),
               p_zip=cbind(psi1, psi2),
               rho = matrix(c(1,.5, .5, 1), nrow=2)
               )

dat = cbind.data.frame(y = dat[[1]],
                       Xt1 = X1[,2],
                       Xt2 = X2[,2],
                       Zt1 = Z1[,2],
                       Zt2 = Z2[,2])

# Generate copula object using copula package
frank.cop <- frankCopula(-2, dim = 2)
gaus.cop <- normalCopula(.5, dim=2)

Xt1 = X1
Xt2 = X2
Zt1 = Z1
Zt2 = Z2

# joint.frank = mvdc(frank.cop,
#                    margins=c("pois", "pois"),
#                    paramMargins = list(list(lambda=as.vector(lam1)),
#                                        list(lambda=as.vector(lam2))
#                    )
# )

# joint.frank = mvdc(frank.cop,
#                    margins=c("zip", "zip"),
#                    paramMargins = list(list(lambda=as.vector(lam1), psi=as.vector(psi1)),
#                                        list(lambda=as.vector(lam2), psi=psi2)
#                    )
# )
# y.frank = rMvdc(n, joint.frank)

# joint.gaus = mvdc(gaus.cop,
#                    margins=c("zinb", "zinb"),
#                    paramMargins = list(list(mu=lam1, size=.5, psi=psi1),
#                                        list(mu=as.vector(lam2), size=.3, psi=psi2)
#                    )
# )

# joint.gaus = mvdc(gaus.cop,
#                    margins=c("zip", "zip"),
#                    paramMargins = list(list(lambda=as.vector(lam1), psi=as.vector(psi1)),
#                                        list(lambda=as.vector(lam2), psi=psi2)
#                    )
#                   )
#
# y.gaus = rMvdc(n, joint.frank)

# y.frank = rMvdc(n, joint.frank)
# y.gaus = rMvdc(n, joint.gaus)

rm(X1, X2, Z1, Z2)

fl1 = y.frank.1 ~ Xt1
fl2 = y.frank.2 ~ Xt2

# z = rzip(n, lambda=exp(2*rnorm(n)), psi=plogis(2*rnorm(n)))
# z = rpois(n, lambda=exp(rnorm(n)))
y1 = rzip(n, lambda=lam1, psi=psi1)
y2 = rzip(n, lambda=lam2, psi=psi2)
# y1 = rpois(n, lambda=lam1)
# y2 = rpois(n, lam2)

dat = cbind.data.frame(y1=y1, y2=y2, Xt1=Xt1[,-1], Xt2=Xt2[,-1], Zt1=Zt1[,-1], Zt2=Zt2[,-1])
# dat = cbind.data.frame(y.frank=y.frank, y.gaus=y.gaus, Xt1=Xt1[,-1], Xt2=Xt2[,-1], Zt1=Zt1[,-1], Zt2=Zt2[,-1], y1, y2)
#
# uni = zic.reg(y.frank.1 ~ Xt1 | Zt1, data=dat, dist="nbinom", optimizer="optim")
# pscl = zeroinfl(y.frank.1~Xt1|Zt1, data=dat, dist="negbin")
#
# cbind(coef(uni), coef(pscl))

test=bizicount(y.1~Xt1|Zt1,
          y.2~Xt2|Zt2,
          margins=c("zip", "zip"),
          cop="gaus",
          start.method="marginal",
          data=dat)
uni1 = coef(zic.reg(y1 ~ Xt1 | Zt1, data = dat))
uni2 = coef(zic.reg(y2 ~ Xt2 | Zt2, data = dat))
# uni1 = coef(glm(y1~Xt1, data=dat, family="poisson"))
# uni2 = coef(glm(y2~Xt2, data=dat, family="poisson"))

# cbind(coef(test), c(b1,b2,NA), c(uni1, uni2, NA))


n = 1000
u = rnorm(n)
x1 = rnorm(n) + u
x2 = rnorm(n) + u
z1 = rnorm(n)
z2 = rnorm(n)
psi1 = plogis(z1*-1)
psi2 = plogis(z2*1.2-1)
lam1 = exp(1.5+x1*2)
lam2 = exp(2+x2*3)

p1 = 1-(psi1 + psi2 - psi1*psi2)
p2 = (1-psi1)*psi2
p3 = (1-psi2)*psi1
p4 = psi1*psi2
y1 = rzip(n, lam1, psi1)
y2 = rzip(n, lam2, psi2)

rbzip.b(n, m0=lam1, m1 = lam2, m2 = exp(rnorm(n)), p1 = p1, p2=p2, p3=p3, p4=p4)

y = rbzip.a(n, m0=lam1, m1 = lam2, m2 = exp(2), p = .25)

zic.reg(y1~x1|z1)

summary(bizicount(y1~x1|z1,
          y2~x2|z2,
          margins = c("zip", "zip"),
          cop="gaus"))

cbind(coef(test), c(b1, b2, g1, g2, NA),
      c(uni1[1:2], uni2[1:2], uni1[3:4], uni2[3:4], NA)
      )

biv = bizicount(y.gaus.1~Xt1+ Zt1|Zt1 + Xt2,
                y.gaus.2~Xt2|Zt2,
                margins = c("zip", "zip"),
                start.method="marginal",
                keep=T,
                data=dat,
                cop="gaus")


gj = gjrm(list(fl1,fl2), BivD="F", Model="B", margins=c("PO", "PO"), data=dat)


cbind(coef(biv), gj$coefficients, biv$se.nid, sqrt(diag(solve(gj$He))))

test=bivcount(y.frank[,1] ~ Xt1 - 1 ,
              y.frank[,2] ~ Xt2 - 1 ,
              margins = c("pois", "pois"),
              start.method="marginal",
              method = c("nlm", "nlminb", "BFGS", "hjn", "Nelder-Mead", "Rvmmin"),
              control = list(maxit=10000),
              keep=T)




joint.frank = mvdc(frank.cop,
                   margins=c("zip", "zip"),
                   paramMargins = list(list(lambda=as.vector(lam1), psi=as.vector(psi1)),
                                       list(lambda=as.vector(lam2), psi=psi2)
                                       )
                   )

joint.gaus = mvdc(gaus.cop,
                  margins=c("zip", "zip"),
                  paramMargins = list(list(lambda=as.vector(lam1), psi=as.vector(psi1)),
                                      list(lambda=as.vector(lam2), psi=psi2)
                                      )
                  )

true.dgp.f = c(b1, b2, g1, g2, disp1=.75, disp2=.5, dep1=.5)
true.dgp.g = c(b1, b2, g1, g2, disp1=.75, disp2=.5, dep1=.5)


y.frank = rMvdc(n, joint.frank)
y.gaus  = rMvdc(n, joint.gaus)

bivcount(y.frank[,1] | y.frank[,2] ~ X - 1,
         margins=c("nbinom", "nbinom"),
         parallel=T, cl.type="FORK")
starts = get.starts("marginal", y.frank, X, Z, c("zinb", "zinb"))
nlm(cop.lik, starts[[2]], X_=X, Z_=Z, y_=y.frank, margins_=c("zinb", "zinb"))

get.cop.starts(y.frank, X, Z, .margins = c("zinb", "zinb"), .cl.type = "PSOCK", .method="stochastic", .parallel=T)
dat   = cbind.data.frame(y=y.frank, x=X[,-1], z=Z[,-1])
dat2  = cbind.data.frame(y=y.gaus, x=X[,-1], z=Z[,-1])
test = bicount(y.1 | y.2 ~ x | z, data=dat, margins=c("zip", "zip"), start.method="naive", cop="frank")
test3 = bicount(y.1 | y.2 ~ x  | z , data=dat2, margins = c("zinb", "zinb"), cop="gaus", start.method = "naive", cl.type="FORK", parallel=T, stepmax=100)
optim(par = c(rep(1,11)), fn=cop.lik, .X=X, Z=Z, y=y.frank,  margins = c("zinb", "zinb"))

x2 = x[, c(1, grep("coords$|pop", colnames(x)))]

test=bicount(y[,"att.ful"] | y[,"att.bok"] ~ x2 + 0 | x2 + 0,
             margins = c("zinb", "zinb"), start.method="naive", cop="gaus", nstarts=6, parallel=T, stepmax=25, keep=T)
test=bicount(X=x2, y=y, Z=x2,
             margins = c("zinb", "zinb"),
             start.method="naive", cop="gaus", nstarts=6, parallel=T, stepmax=25, keep=T)
test2=bicount(y[,"att.ful"] | y[,"att.bok"] ~ x2 + 0 | x2 + 0,
              start.method="naive", cop="gaus", nstarts=22, parallel=T, stepmax=5)

load("./data/data_john.RData")

pscl = zeroinfl(y[,1] ~ x + 0 | x+0, dist="negbin")
zic = zic.reg(y[,1] ~ x + 0, ~ x+ 0, dist="negbin", optimizer = "optim", method="Nelder-Mead")

dat2 = dat[dat$xcoord != 3.25 | dat$ycoord != 6.75, ]
dat2 = dat[,-which(names(dat) %in% c("xcoord", "ycoord"))]
dat3 = dat2[,c("att.bok", "att.ful", "pop_gpw_sum_2010", "xcoords", "ycoords", "mountains_mean")]

dat3 = dat2[, -grep("coords$", names(dat2))]
start = Sys.time()
emp.f = tryCatch(bicount(att.bok |
                  att.ful ~ .  | . ,
                data=dat3, margins=c("zip", "zinb"), cop="frank",
                keep=T, parallel=T),
                error = function(e) conditionMessage(e))
Sys.time() - start

x2 = x[, c( "pop_gpw_sum_2010", "xcoords", "ycoords", "mountains_mean")]
x2 = cbind(int=1, x2)


starts.n = starts.naive(y, x2, x2, c("zinb", "zinb"))
starts.m = starts.marginal(y, x2, x2, c("zinb", "zinb"))

test=nlm(p=starts.m, f = cop.lik, X_=x2, Z_=x2, y_=y, margins=c("zinb", "zinb"), iterlim=500, gradtol=1e-2)
test2=nlm(p=starts.n, f = cop.lik, X_=x2, Z_=x2, y_=y, margins=c("zinb", "zinb"), iterlim=500, gradtol=1e-2)

nlm(p=test2$estimate, f = cop.lik, X_=x2, Z_=x2, y_=y, margins=c("zinb", "zinb"), iterlim=500)

start = Sys.time()
emp.f2 = tryCatch(bicount(att.ful |
                           att.bok ~ .  | . ,
                         data=dat3, margins=c("zip", "zip"), cop="frank",
                         keep=T),
                 error = function(e) conditionMessage(e))

Sys.time() - start

emp.f3 = tryCatch(bicount(att.ful |
                           att.bok ~ .  | . ,
                         data=dat3, margins=c("zip", "zip"), cop="frank",
                         keep=T, parallel=F, starts=rep(1, 21)),
                 error = function(e) conditionMessage(e))



emp.f2 = bicount(att.ful |
                   att.bok ~ .  | . ,
                 data=dat2, margins=c("zip", "zip"), cop="frank", keep=T, start.method = "marginal", nstarts=6, parallel=T)
emp.f3 = bicount(att.ful |
                  att.bok ~ .  | . ,
                data=dat2, margins=c("zip", "zip"), cop="frank", keep=T, start.method = "stochastic", nstarts=6, parallel=T)

test.frank = nlm(f=cop.lik, p=rep(1, 11), X=X, y=y.frank, Z=Z, cop="frank", margins=c("zinb", "zinb"), iterlim=5000, pmf.min = 1e-7, hessian=T, offset.ct=0, offset.zi=0)

nloptr(eval_f=cop.lik, x0=rep(1, 11), X=X, y=y.frank, Z=Z,
    cop="frank", margins=c("zinb", "zinb"), pmf.min = 1e-7, bounds=T, frech.min=1e-7, opts=list(algorithm="NLOPT_LN_SBPLX", maxeval=5000 ))

start1 = glm.fit( X, y.gaus[,1], family=poisson())$coefficients
start2 = glm.fit(X, y.gaus[,2], family=poisson())$coefficients
start3 = glm.fit( Z, as.integer(y.gaus[,1]==0), family=binomial())$coefficients
start4 = glm.fit( Z, as.integer(y.gaus[,2]==0), family=binomial())$coefficients
start.disp = c(.01, .01)
start.dep = .01

starts = c(start1, start2, start3, start4, start.disp, start.dep)

# est1 = c(231.039, 335.8784, 657.6105, 795.7263)
# est2 = c(-794.2175, 23.99494, -783.7254, 15.34981)
# disp = exp(c(131.4057, 59.22961))
# dep = tanh(2004.738)
# rm(est1, est2, disp, dep)


test.gaus = nlm(f=cop.lik, p=rep(1,11), X=X, y=y.gaus, Z=Z, cop="gaus",
                margins=c("zinb", "zinb"), iterlim=5000, pmf.min = 1e-7, hessian=T, stepmax = 25)


cbind(test.gaus$estimate, test.frank$estimate, true.dgp)
#=================
# zinb pois

joint.frank = mvdc(frank.cop,
                   margins=c("zinb", "pois"),
                   paramMargins = list(list(mu=as.vector(lam1), size = .75, psi=as.vector(psi1)),
                                       list(lambda=as.vector(lam2))
                                       )
                   )
joint.gaus = mvdc(gaus.cop,
                  margins=c("zinb", "pois"),
                  paramMargins = list(list(mu=as.vector(lam1), size = .75, psi=as.vector(psi1)),
                                      list(lambda=as.vector(lam2))
                                      )
                  )

y.frank = rMvdc(n, joint.frank)
y.gaus  = rMvdc(n, joint.gaus)

start1 = glm.fit( X, y.gaus[,1], family=poisson())$coefficients
start2 = glm.fit(X, y.gaus[,2], family=poisson())$coefficients
start3 = glm.fit( Z, as.integer(y.gaus[,1]==0), family=binomial())$coefficients

start.disp = c(.05)
start.dep = 0

starts = c(start1,start2,start3,start.disp, start.dep)
test.frank = nlm(f=cop.lik, p=rep(1, 8), X=X, y=y.frank, Z=Z, cop="frank", margins=c("zinb", "pois"), iterlim=5000, pmf.min = 1e-7, hessian=T)

test.gaus = nlm(f=cop.lik, p=rep(1,8), X=X, y=y.gaus, Z=Z, cop="gaus",
                margins=c("zinb", "pois"), iterlim=5000, pmf.min = 1e-7, hessian=T, stepmax=25)


#================
# zip nbinom

joint.frank = mvdc(frank.cop,
                   margins=c("zip", "nbinom"),
                   paramMargins = list(list(lambda=as.vector(lam1), psi=as.vector(psi1)),
                                       list(mu=as.vector(lam2), size=.75)
                   )
)

joint.gaus = mvdc(gaus.cop,
                  margins=c("zip", "nbinom"),
                  paramMargins = list(list(lambda=as.vector(lam1), psi=as.vector(psi1)),
                                      list(mu=as.vector(lam2), size = .75)
                  )
)

y.frank = rMvdc(n, joint.frank)
y.gaus  = rMvdc(n, joint.gaus)

start1 = glm.fit( X, y.gaus[,1], family=poisson())$coefficients
start2 = glm.fit(X, y.gaus[,2], family=poisson())$coefficients
start3 = glm.fit( Z, as.integer(y.gaus[,1]==0), family=binomial())$coefficients

start.disp = 0.01
start.dep = 0

starts = c(start1,start2,start3, start.disp, start.dep)
test.frank = nlm(f=cop.lik, p=rep(1, 8), X=X, y=y.frank, Z=Z, cop="frank", margins=c("zip", "nbinom"), iterlim=5000, pmf.min = 1e-7, hessian=T)

test.gaus = nlm(f=cop.lik, p=rep(1, 8), X=X, y=y.gaus, Z=Z, cop="gaus",
                margins=c("zip", "nbinom"), iterlim=5000, pmf.min = 1e-7, hessian=T, stepmax=25)



dat=rzip(100, runif(100, 0, 5), psi=runif(100))
dat=rzip(100, runif(100, -1, 5), psi=runif(100))

mfn = function(cd=NULL){
  assert(!is.null(cd), "an error", type="message", appendLF=F)
}
mfn()
for(i in 1:5){mfn();Sys.sleep(.1)}
cond = structure(list(message="test", call=NULL), class=c("warning", "condition"))
for(i in 1:5){warning(cond);Sys.sleep(.5)}

rzinb(10, mu=10, size = -1, psi=.5)
check = function() parent.frame(1)

my.call = "A call"

test = tryCatch(warning("A warning", immediate.=T),
         warning = function(w) w
)

for(i in 1:20){

  tryCatch(warning("A warning", immediate.=F),
           warning = function(w){w$call <- my.call; warning(w)}
           )
  Sys.sleep(.1)

}


microbenchmark({rzip(100, runif(100, 0, 5), psi=runif(100))}, unit="s")


n=1000
X=rnorm(n)
Z = runif(n)
X = cbind(1, X)
Z = cbind(1, Z)
b = c(.3, 1.5)
g = c(1, -3)

psi = plogis(Z%*%g)
lam = exp(X%*%b)

y = rzinb(n, mu=lam, psi=psi, size=.6)
pois.grad(c(b,g), X=X, z=Z, y=y)
glm.fit(X, y, family=poisson(link="sqrt"))

dat= cbind.data.frame(X=X[,2], Z=Z[,2], y=y)
test = zic.reg(y~X + offset(offset.ct) |Z , dist="nbinom", data=dat)
test2 = zeroinfl(y~X+0|Z+0, dist="negbin", link="cloglog")
nlm(f=zic.ll, p=rep(1,5), X=X, z=Z, y=y, dist="nbinom", grad=T)

library(microbenchmark)

zic.ll = cmpfun(zic.ll)
nbinom.grad = cmpfun(nbinom.grad)
pois.grad = cmpfun(pois.grad)
dzinb = cmpfun(dzinb)


post=microbenchmark(setup = {n=1000;
                        X=rnorm(n);
                        Z = runif(n);
                        X = cbind(1, X);
                        Z = cbind(1, Z);
                        b = c(-2, 1.5);
                        g = c(1, -3);
                        psi = plogis(Z%*%g);
                        lam = exp(X%*%b);
                        y = rzinb(n, mu=lam, psi=psi, size=.64)
                        },
               op.bf = zic.reg(y~X+0|Z+0, dist="nbinom", optimizer="optim"),
               op.n = zic.reg(y~X+0|Z+0, dist="nbinom", optimizer="optim", method="Nelder-Mead"),
               nlm = zeroinfl(y~X+0|Z+0, dist="negbin", optimizer="nlm")
)

psi = plogis(Z%*%g)
lam = exp(X%*%b)

y = rzinb(n, lambda=lam, psi=psi, size=.64)


test=nbinom.grad(rep(1,5), X=x, z=Z, y=y, offset.ct=0, offset.zi=0, weights=1)
gradNegBin(rep(1,5))

optim(par=rep(1,4), fn=zic.ll, X=X, y=y, z=Z, dist="pois", method="BFGS", gr=pois.grad)
nlm(rep(1,5), f=zic.ll, X=X, y=y, z=Z, dist="nbinom",  grad=T )

# discResidFn = function(data){
#   mn =  pzinb(data-1, mu=lam2, psi=psi2, size=.7)
#   mx =  pzinb(data, mu=lam2, psi=psi2, size=.7)
#
#   out = ifelse(mn==mx, mn, runif(length(data), mn, mx))
#
#   return(out)
#   out = vector(length=length(mn))
#   out[mn==mx] = mn
#   out[mn!=mx] = runif(1, mn, mx)
#   return(out)
#
#   return(list(mn, mx))
#   if(mn==mx)
#     return(mn)
#   else
#     return(
#       runif(1, mn, mx)
#     )
# }
