# zi poisson pmf ===============================================================
na_warn = "NaNs produced; perhaps `lambda` negative, or `psi` not in [0, 1]."
na_warn_nb = "NaNs produced; perhaps `mu` or `size` negative, or `psi` or `prob` not in [0, 1]?"

#' @name zip
#' @title The zero-inflated Poisson (ZIP) distribution
#'
#' @description These functions are used to evaluate the zero-inflated Poisson
#'   distribution's probability mass function (PMF), cumulative distribution
#'   function (CDF), and quantile function (inverse CDF), as well as generate
#'   random realizations from the ZIP distribution.
#'
#' @details The zero inflated Poisson distribution is a mixture of a Poisson
#'    and a degenerate point-mass at 0. It has the form
#'    \deqn{\psi + (1-\psi)(\lambda^x e^-\lambda)/x!}, with mean
#'    \eqn{(1-\psi)\lambda}. Thus, the parameter `lambda` above is the mean of
#'    the Poisson distribution that forms part of the zero-inflated
#'    distribution, *not* the mean of the ZIP distribution.
#'
#'    `recycle` -- If `FALSE` (default), all arguments must have identical
#'    length, there can be two unique lengths for the arguments, provided that
#'    one of those lengths is 1. For example, `lambda = c(1,2,3)` and `psi=.5`
#'    is acceptable because there are two unique lengths, and one of them is
#'    length 1. However, `lambda=c(1,2,3)` and `psi=c(.5,.2)` would fail, as
#'    there are two distinct lengths, none of which is 1. If `TRUE,` no
#'    additional checks (beyond those in base R's functions) are made to ensure
#'    that the argument vectors have the same length.
#'
#' @param x,q Vector of quantiles at which to evaluate the PMF and CDF,
#'  respectively. Should be non-negative integers.
#'
#' @param p Vector of probabilities at which to evaluate the quantile function.
#'
#' @param n Number of realizations from the distribution to generate
#'
#' @param lambda Vector of means for the count portion of the zero-inflated
#'    Poisson distribution. Should be non-negative. NOTE: This is *not* the mean
#'    of the zero-inflated Poisson distribution; it is the mean of the Poisson
#'    component of the mixture distribution. See 'Details.'
#'
#' @param psi Vector of zero-inflation probabilities.
#'
#' @param log,log.p Logical indicating whether probabilities should be returned
#'  on log scale (for `dzip` and `pzip`), or are supplied on log-scale (for `qzip`).
#'
#' @param lower.tail Logical indicating whether probabilities should be
#'    \eqn{Pr(X \le x)} or \eqn{Pr(X > x)}
#'
#' @param recycle Logical indicating whether to permit arbitrary recycling of
#'    arguments with unequal length. See 'Details' and 'Examples.'
#'
#' @return `dzip` returns the mass function evaluated at `x`,
#'    `pzip` returns the CDF evaluated at `q`, `qzip` returns the quantile
#'    function evaluated at `p`, and `rzip` returns random variates with the
#'    specified parameters.
#'
#' @example inst/examples/zip_dist_ex.R
#'
#' @author John Niehaus
#'
#' @references Lambert, Diane. "Zero-inflated Poisson regression, with an
#'   application to defects in manufacturing." Technometrics 34.1 (1992): 1-14.
#'
#' @export
dzip = function(x, lambda, psi, log = FALSE, recycle = FALSE){

     if(!env_has(e.check, "zi.dens.checks"))
          check_dist_args(recycle=recycle)

     call = match.call()

     l = log
     rm(log)

     n = max(length(x), length(lambda), length(psi))

     if( all(lambda < 0) || !any(vprob(psi))){
          warning(na_warn)
          return(rep(NaN, n))
     }

     new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
     list2env(new.args, envir=environment())

     # note that doing logs this way is computationally superior, as it uses the internal logarithm
     # for the dpois function and log1p for small psi/1-psi
     # computing log format first, then exponentiating if necessary also avoids issues with log(0)
     # because sometimes 0 occurs due to underflow
     x = rep(x, length.out=n)
     x0 = x == 0
     pmf = alter.cond(log1p(-psi) + dpois(x = x, lambda = lambda, log = T), suppress = T)
     pmf[x0] = suppressWarnings(log(psi + exp(pmf))[x0])

     pmf[psi < 0 | psi > 1] = NaN

     if(anyNA(pmf))
          warning(na_warn)

     if(!l)
          pmf = exp(pmf)

     return(pmf)

}

# zi poisson rng ===============================================================
#' @rdname zip
#' @export
rzip = function(n, lambda, psi, recycle = FALSE){

     check_dist_args(recycle=recycle)

     call = match.call()

     if( all(lambda < 0) || !any(vprob(psi))){
          warning(na_warn)
          return(rep(NaN, n))
     }

     new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
     list2env(new.args, envir=environment())

     ind = suppressWarnings(rbinom(n, 1, psi))
     x = suppressWarnings(rpois(n, lambda=lambda))

     x[ind == 1 & !is.na(x)] = 0
     x[is.na(ind) | is.na(x)] = NaN
     if(anyNA((x)))
          warning(na_warn)

     return(x)
}



# zi poisson cdf ===============================================================
#' @rdname zip
#' @export
pzip = function(q, lambda, psi, lower.tail = TRUE, log.p = FALSE, recycle = FALSE){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(recycle=recycle)
  }

  call = match.call()

  n = max(length(q), length(lambda), length(psi))
  if( all(lambda < 0) || !any(vprob(psi))){
     warning(na_warn)
     return(rep(NaN, n))
  }

  new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  cdf = suppressWarnings(psi + (1-psi)*ppois(q, lambda))



  cdf[q < 0] = 0
  cdf[psi > 1 | psi < 0 | lambda < 0] = NaN

  if(anyNA(cdf))
       warning(na_warn)

  if(!lower.tail)
    cdf = 1 - cdf

  if(log.p)
    cdf = log(cdf)

  return(cdf)

}


# zi poisson quantile =========================================================
#' @rdname zip
#' @export
qzip = function(p, lambda, psi, lower.tail = TRUE, log.p = FALSE, recycle = FALSE){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(recycle=recycle)
  }

  call = match.call()

  n = max(length(p), length(lambda), length(psi))

  if( all(lambda < 0) || !any(vprob(psi)) || !any(vprob(p)) ){
    warning(na_warn)
    return(rep(NaN, n))
  }

  new.args = lapply(list(psi=psi, lambda=lambda), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  if(!lower.tail)
    p = 1 - p

  if(log.p)
    p = exp(p)

  p2 = (p - psi)/(1 - psi)

  quant = suppressWarnings( qpois(p2, lambda = lambda) )

  quant[p < psi] = 0
  quant[psi < 0 | psi > 1 | p < 0 | p > 1 | lambda < 0] = NaN

  if(anyNA(quant)){
    warning(na_warn)
  }

  return(quant)
}



### PMF, CDF, Quantile, Random generation for zero inflated negbin

# zi nbin pmf =================================================================

#' @name zinb
#' @title The zero-inflated negative binomial (ZINB) distribution
#'
#' @description These functions are used to evaluate the zero-inflated negative binomial
#'   distribution's probability mass function (PMF), cumulative distribution
#'   function (CDF), and quantile function (inverse CDF), as well as generate
#'   random realizations from the ZINB distribution.
#'
#' @param x,q Vector of quantiles at which to evaluate the PMF and CDF,
#'  respectively. Should be non-negative integers.
#'
#' @param p Vector of probabilities at which to evaluate the quantile function.
#'
#' @param n Number of realizations to generate from the distribution
#'
#' @param mu Vector of means for the count portion of the zero-inflated negative
#'   binomial distribution. Only one of `mu` or `prob` should be specified, not
#'   both. Should be non-negative. NOTE: This is *not* the mean of the ZINB
#'   distribution; it is the mean of the NB component of the mixture
#'   distribution. See \code{\link[stats:dnbinom]{nbinom}}.
#'
#' @param size The inverse dispersion parameter, or number of successful trials,
#'   both for the negative binomial portion of the ZINB mixture distribution.
#'   See \code{\link[stats:dnbinom]{nbinom}}.
#'
#' @param prob The probability of success on each trial in the negative binomial portion of the mixture distribution. Only one of `mu` or
#'   `prob` should be specified, not both. See \code{\link[stats:dnbinom]{nbinom}}.
#'
#' @param psi Vector of zero-inflation probabilities.
#'
#' @param log,log.p Logical indicating whether probabilities should be returned
#'  on log scale (for `dzip` and `pzip`), or are supplied on log-scale (for `qzip`).
#'
#' @param lower.tail Logical indicating whether probabilities should be
#'    \eqn{Pr(X \le x)} or \eqn{Pr(X > x)}
#'
#' @param recycle Logical indicating whether to permit arbitrary recycling of
#'    arguments with unequal length. See 'Details' and 'Examples.'
#'
#' @return `dzinb` returns the mass function evaluated at `x`,
#'    `pzinb` returns the CDF evaluated at `q`, `qzinb` returns the quantile
#'    function evaluated at `p`, and `rzinb` returns random realizations with the
#'    specified parameters.
#'
#' @author John Niehaus
#'
#' @references Lambert, Diane. "Zero-inflated Poisson regression, with an
#'   application to defects in manufacturing." Technometrics 34.1 (1992): 1-14.
#'
#' @example /inst/examples/zinb_dist_ex.R
#'
#' @export
dzinb = function(x, size, psi, mu = NULL, prob = NULL, lower.tail = TRUE, log = FALSE, recycle = FALSE){

     if(!env_has(e.check, "zi.dens.checks")){
          check_dist_args(negbin=T, recycle=recycle)
     }

     call = match.call()

     l = log
     rm(log)

     if(is.null(mu))
          mu = size*(1-prob)/prob

     n = max(length(x), length(size), length(psi), length(prob), length(mu))
     if (all(size < 0) || !any(vprob(psi)) || (!is.null(mu) && all(mu < 0)) || !any(vprob(prob))) {
          warning(na_warn_nb)
          return(rep(NaN, n))
     }

     x = rep(x, length.out = n)
     new.args = lapply(list(size=size, psi=psi, mu=mu, prob=prob), rep_len.custom, n=n)
     list2env(new.args, envir=environment())

     x0 = x == 0
     pmf =   alter.cond(log1p(-psi) + dnbinom(x=x, size=size, mu=mu, log=T), suppress=T)
     pmf[x0] = suppressWarnings( log(psi + exp(pmf))[x0] )

     pmf[psi < 0 | psi > 1 | prob > 1 | prob <= 0] = NaN

     if(!l)
          pmf = exp(pmf)

     if(anyNA(pmf))
          warning(na_warn_nb)

     return(pmf)

}



# zi nbin cdf =================================================================
#' @rdname zinb
#' @export
pzinb = function(q, size, psi,  mu = NULL, prob = NULL, lower.tail = TRUE, log.p = FALSE, recycle = FALSE){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(negbin=T, recycle=recycle)
  }

  if(is.null(mu))
    mu = size*(1-prob)/prob

  n = max(length(q), length(size), length(psi), length(prob), length(mu))
  if (all(size < 0) || !any(vprob(psi)) || (!is.null(mu) && all(mu < 0)) || !any(vprob(prob, eq = 'right'))) {
       warning(na_warn_nb)
       return(rep(NaN, n))
  }


  new.args = lapply(list(size=size, psi=psi, mu=mu, prob=prob), rep_len.custom, n=n)
  list2env(new.args, envir=environment())

  call = match.call()

  cdf = alter.cond(
      psi + (1-psi)*pnbinom(q=q, size=size, mu=mu, log.p=F), suppress=T
  )

  cdf[q<0] = 0
  cdf[mu < 0 | psi < 0 | psi > 1 | prob > 1 | prob <= 0 | size < 0] = NaN

  if(!lower.tail)
    cdf = 1 - cdf

  if(log.p)
    cdf = log(cdf)

  if(anyNA(cdf))
    warning(na_warn_nb)

  return(cdf)

}


# zi nbin quantile =============================================================
#' @rdname zinb
#' @export
qzinb = function(p, size, psi, mu = NULL, prob = NULL, lower.tail = TRUE, log.p = FALSE, recycle = FALSE){

  if(!env_has(e.check, "zi.dens.checks")){
    check_dist_args(negbin=T, recycle=recycle)
  }

  if(is.null(mu))
    mu = size*(1-prob)/prob

  n = max(length(p), length(size), length(psi), length(prob), length(mu))
  if (all(size < 0) || !any(vprob(psi)) || (!is.null(mu) && all(mu < 0)) || !any(vprob(prob, eq = 'right')) || !any(vprob(p))) {
       warning(na_warn_nb)
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
  quant[mu < 0 | p < 0 | p > 1 | psi < 0 | psi > 1 | prob > 1 | prob <= 0  | size < 0] = NaN

  if(anyNA(quant))
    warning(na_warn_nb)

  return(quant)
}


# zi nbin rng =================================================================
#' @rdname zinb
#' @export
rzinb = function(n, size, psi, mu = NULL, prob = NULL, recycle = FALSE){

     check_dist_args(negbin=T, recycle=recycle)
     call = match.call()

     if (all(size < 0) || !any(vprob(psi)) || (!is.null(mu) && all(mu < 0)) || !any(vprob(prob, eq = 'right'))) {
          cat("here")
          warning(na_warn_nb)
          return(rep(NaN, n))
     }


     if(is.null(mu))
          mu =  size*(1-prob)/prob

     new.args = lapply(list(size=size, psi=psi, mu=mu, prob=prob), rep_len.custom, n=n)
     list2env(new.args, envir=environment())

     ind = suppressWarnings(rbinom(n, 1, psi) == 1)
     x = suppressWarnings(rnbinom(n, size = size, mu = mu))

     x[ind == 1 & !is.na(x)] = 0
     x[is.na(ind) | is.na(x)] = NaN

     if(anyNA(x))
          warning(na_warn_nb)

     return(x)
}



