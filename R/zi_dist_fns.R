# zi poisson pmf ===============================================================
dzip = function(x, lambda, psi, log=F, recycle=F){

  if(!env_has(e.check, "zi.dens.checks"))
    check_dist_args(recycle=recycle)

  call = match.call()

  l = log
  rm(log)

  n = max(length(x), length(lambda), length(psi))

    if((length(lambda)==1 && lambda < 0) ||
     (length(psi) ==1 && (psi < 0 || psi > 1))){
    warning("NaNs produced; perhaps `lambda` negative, or `psi` not on unit-interval?")
    return(rep(NaN, n))
  }

  new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  # note that doing logs this way is computationally superior, as it uses the internal logarithm
  # for the dpois function and log1p for small psi/1-psi
  # computing log format first, then exponentiating if necessary also avoids issues with log(0)
  # because sometimes 0 occurs due to underflow
  x = rep(x, length.out=n)
  pmf = alter.cond(
      ifelse(
        x==0,
        log(psi +  exp( log1p(-psi) + dpois(x=x, lambda=lambda, log=T) )),
        log1p(-psi) + dpois(x=x, lambda=lambda, log=T)),
      new.message = "NaNs produced; perhaps `lambda` negative, or `psi` not on unit-interval?"
  )

  pmf[lambda < 0 | psi < 0 | psi > 1] = NaN

  if(!l)
    pmf = exp(pmf)

  return(pmf)

}


# zi poisson cdf ===============================================================
pzip = function(q, lambda, psi, lower.tail=T, log.p=F, recycle=F){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(recycle=recycle)
  }

  call = match.call()

  n = max(length(q), length(lambda), length(psi))
  if((length(lambda)==1 && lambda < 0) ||
     (length(psi) ==1 && (psi < 0 || psi > 1))){
    warning("NaNs produced; perhaps `lambda` negative, or `psi` not on unit-interval?")
    return(rep(NaN, n))
  }

  new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  cdf = alter.cond(
      psi + (1-psi)*ppois(q, lambda),
      new.message = "NaNs produced; perhaps `lambda` negative, or `psi` not on unit-interval?"
  )

  cdf[lambda < 0 | psi > 1 | psi <0] = NaN
  cdf[q < 0] = 0

  if(!lower.tail)
    cdf = 1 - cdf

  if(log.p)
    cdf = log(cdf)

  return(cdf)

}


# zi poisson quantile =========================================================
qzip = function(p, lambda, psi, lower.tail=T, log.p=F, recycle=F){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(recycle=recycle)
  }

  call = match.call()

  n = max(length(p), length(lambda), length(psi))

  if(all(lambda < 0) || all(psi < 0 | psi > 1) || all(p < 0)){
    warning("NaNs produced; perhaps `lambda` negative, or `p`, `psi` not on unit-interval?")
    return(rep(NaN, n))
  }

  new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  if(!lower.tail)
    p = 1 - p

  if(log.p)
    p = exp(p)

  p2 = (p - psi)/(1 - psi)

  quant= alter.cond(
    qpois(p2, lambda = lambda),
    suppress = T
  )

  quant[p < psi] = 0
  quant[lambda < 0 | psi< 0 | psi > 1 | p< 0 | p>1] = NaN

  if(any(is.na(quant))){
    warning("NaNs produced; perhaps `lambda` negative, or `p`, psi` not on unit-interval?")
  }

  return(quant)
}



# zi poisson rng ===============================================================
rzip = function(n, lambda, psi, recycle=F){

  check_dist_args(recycle=recycle)

  call = match.call()

  if((length(lambda)==1 && lambda < 0) ||
     (length(psi) ==1 && (psi < 0 || psi > 1))){
    warning("NaNs produced; perhaps `lambda` negative, or `psi` not on unit-interval?")
    return(rep(NaN, n))
  }

  new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  x = alter.cond(
    ifelse(rbinom(n,1,psi) == 1, 0, rpois(n, lambda=lambda))
  )

  x[psi <0 | psi> 1 | lambda < 0] = NaN

  return(x)
}




### PMF, CDF, Quantile, Random generation for zero inflated negbin

# zi nbin pmf =================================================================
dzinb = function(x, size, psi,  prob=NULL, mu=NULL, log=F, recycle=F){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(negbin=T, recycle=recycle)
  }

  call = match.call()

  l = log
  rm(log)

  if(is.null(mu))
    mu = size*(1-prob)/prob

  n = max(length(x), length(size), length(psi), length(prob), length(mu))
  if((length(size) == 1 && size <= 0) ||
     (length(psi) == 1 && (psi > 1 | psi < 0)) ||
     (length(mu) == 1 && mu < 0) ||
     length(prob) == 1 && (prob > 1 | prob <= 0)
  ){
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi` or `prob` not on unit-interval?")
    return(rep(NaN, n))
  }

  x = rep(x, length.out = n)
  new.args = lapply(list(size=size, psi=psi, mu=mu, prob=prob), rep_len.custom, n=n)
  list2env(new.args, envir=environment())


  pmf = alter.cond(
        ifelse(
          x==0,
          log(psi +  exp( log1p(-psi) + dnbinom(x=x, size=size, mu=mu, log=T) )),
          log1p(-psi) + dnbinom(x=x, size=size, mu=mu, log=T)
        )
  )

  if(!l)
    pmf = exp(pmf)

  pmf[mu<0|psi < 0 | psi > 1 | prob > 1 | prob <= 0 | size <= 0] = NaN

  if(anyNA(pmf))
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi`, `prob` not on unit-interval?")

  return(pmf)

}



# zi nbin cdf =================================================================
pzinb = function(q, size, psi, prob=NULL, mu=NULL,  lower.tail=T, log.p=F, recycle=F){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(negbin=T, recycle=recycle)
  }

  if(is.null(mu))
    mu = size*(1-prob)/prob

  n = max(length(q), length(size), length(psi), length(prob), length(mu))
  if((length(size) == 1 && size <= 0) ||
     (length(psi) == 1 && (psi > 1 | psi < 0)) ||
     (length(mu) == 1 && mu < 0) ||
     length(prob) == 1 && (prob > 1 | prob <= 0)
  ){
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi` or `prob` not on unit-interval?")
    return(rep(NaN, n))
  }

  new.args = lapply(list(size=size, psi=psi, mu=mu, prob=prob), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  call = match.call()

  cdf = alter.cond(
      psi + (1-psi)*pnbinom(q=q, size=size, mu=mu, log.p=F)
  )

  cdf[q<0] = 0
  cdf[mu < 0 | psi < 0 | psi > 1 | prob > 1 | prob <= 0 | size <= 0] = NaN

  if(!lower.tail)
    cdf = 1 - cdf

  if(log.p)
    cdf = log(cdf)

  if(anyNA(cdf))
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi`, `prob` not on unit-interval?")

  return(cdf)

}


# zi nbin quantile =============================================================
qzinb = function(p, size, psi, prob=NULL, mu=NULL, lower.tail=T, log.p=F, recycle=F){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(negbin=T, recycle=recycle)
  }

  if(is.null(mu))
    mu = size*(1-prob)/prob

  n = max(length(p), length(size), length(psi), length(prob), length(mu))
  if(all(size <= 0) ||
     all(psi > 1 | psi < 0) ||
     all(mu < 0) ||
     (all(prob > 1 | prob <= 0) && !is.null(prob))
  ){
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi` or `prob` not on unit-interval?")
    return(rep(NaN, n))
  }

  new.args = lapply(list(size=size, psi=psi, mu=mu, prob=prob), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  call = match.call()

  if(log.p)
    p = exp(p)

  if(!lower.tail)
    p = 1 - p

  p2 = (p - psi)/(1 - psi)

  quant = alter.cond(
    qnbinom(p2, size = size, mu = mu),
    suppress=T
  )

  quant[p < psi] = 0
  quant[mu<0| p<0 | p > 1 | psi < 0 | psi > 1 | prob > 1 | prob <= 0 | size <= 0] = NaN

  if(anyNA(quant))
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi`, `prob`, `q` not on unit-interval?")

  return(quant)
}


# zi nbin rng =================================================================
rzinb = function(n, size, psi, mu=NULL, prob=NULL, recycle=F){

  check_dist_args(negbin=T, recycle=recycle)
  call = match.call()

  if((length(size) == 1 && size <= 0) ||
     (length(psi) == 1 && (psi > 1 | psi < 0)) ||
     (length(mu) == 1 && mu < 0) ||
     length(prob) == 1 && (prob > 1 | prob <= 0)
  ){
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi` or `prob` not on unit-interval?")
    return(rep(NaN, n))
  }

  if(is.null(mu))
    mu =  size*(1-prob)/prob

  new.args = lapply(list(size=size, psi=psi, mu=mu, prob=prob), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  x = alter.cond(
      ifelse(rbinom(n,1,psi) == 1, 0, rnbinom(n, size=size, mu=mu))
  )

  x[mu < 0 | psi < 0 | psi > 1 | prob > 1 | prob <= 0 | size <= 0] = NaN
  if(anyNA(x))
    warning("NaNs produced; perhaps `mu` or `size` negative, or `psi` or `prob` not on unit-interval?")

  return(x)
}


