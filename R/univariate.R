#' @importFrom rlang env_has
#' @importFrom rlang env_poke
#' @import Formula

#univariate regression function
zic.reg = function(fmla = NULL,
                   data,
                   subset,
                   na.action,
                   weights = rep(1, length(y)),
                   X = NULL,
                   z = NULL,
                   y = NULL,
                   offset.ct = NULL,
                   offset.zi = NULL,
                   dist = "pois",
                   link.ct = "log",
                   link.zi = "logit",
                   starts = NULL,
                   optimizer = "nlm",
                   warn.parent = T,
                   keep = F,
                   ...) {

  if (!env_has(e.check, "univ.check"))
    check_uni_args()

  # Univariate Likelihood
  zic.ll = function(parms,
                    grad = F) {
    kx = ncol(X)
    kz = ncol(z)
    #setup parameters based on distribution
    if (dist == "pois") {
      ct = head(parms, kx)
      zi = tail(parms, kz)
    } else {
      ct = head(parms, kx)
      zi = parms[(kx + 1):(kx + kz)]
      disp = exp(tail(parms, 1))
    }

    #get linear predictor, probability, mean
    eta.ct = X %*% ct + offset.ct
    eta.zi = z %*% zi + offset.zi
    psi = make.link(link.zi)$linkinv(eta.zi)
    lam = make.link(link.ct)$linkinv(eta.ct)

    ll = withCallingHandlers(switch(
      dist,
      "pois" = dzip(
        y,
        lambda = lam,
        psi = psi,
        log = T
      ),
      "nbinom"  = dzinb(
        y,
        size = disp,
        psi = psi,
        mu = lam,
        log = T
      )
    ),
    warning = ll.suppress)

    ll =  sum(ll * weights)

    if (grad)
      attr(ll, "gradient") = get(paste0(dist, ".grad"))(parms = parms)

    return(-ll)
  }

  # Univariate zip gradient
  pois.grad = function(parms,
                       grad = F) {

    kx = ncol(X)
    kz = ncol(z)

    ct = head(parms, kx)
    zi = tail(parms, kz)

    eta.ct = (X %*% ct) + offset.ct
    eta.zi = (z %*% zi) + offset.zi

    linklist.ct = make.link(link.ct)
    linklist.zi = make.link(link.zi)

    lam.ct = linklist.ct$linkinv(eta.ct)
    deta.ct = linklist.ct$mu.eta(eta.ct)

    psi = linklist.zi$linkinv(eta.zi)
    deta.zi = linklist.zi$mu.eta(eta.zi)


    denom.y0 = psi + (1 - psi) * exp(-lam.ct)

    d.zi = ifelse(y == 0,
                  deta.zi * (1 - exp(-lam.ct)) / denom.y0,
                  -deta.zi / (1 - psi))

    d.ct = ifelse(y == 0,
                  -(1 - psi) * exp(-lam.ct) * deta.ct / denom.y0,
                  deta.ct * (y / lam.ct - 1))

    d.zi = z * d.zi * weights
    d.ct = X * d.ct * weights


    return(-unname(colSums(cbind(d.ct, d.zi))))
  }


  #univariate zinb gradient
  nbinom.grad = function(parms,
                         grad = F) {
    kx = ncol(X)
    kz = ncol(z)
    ct = head(parms, kx)
    zi = parms[(kx + 1):(kx + kz)]
    theta = exp(tail(parms, 1))

    eta.ct = (X %*% ct) + offset.ct
    eta.zi = (z %*% zi) + offset.zi

    linklist.ct = make.link(link.ct)
    linklist.zi = make.link(link.zi)

    lam.ct = linklist.ct$linkinv(eta.ct)
    deta.ct = linklist.ct$mu.eta(eta.ct)

    psi = linklist.zi$linkinv(eta.zi)
    deta.zi = linklist.zi$mu.eta(eta.zi)

    denom.y0 = dzinb(y, size = theta, psi = psi, mu = lam.ct) ^ (-1)
    dnbinom.zi = dnbinom(y, size = theta, mu = lam.ct)

    d.zi = ifelse(y == 0,
                  denom.y0 * (deta.zi - deta.zi * dnbinom.zi),
                  -deta.zi / (1 - psi))

    d.ct = ifelse(
      y == 0,
      -denom.y0 * deta.ct * (1 - psi) * (theta / (lam.ct + theta)) ^ (theta + 1),
      deta.ct * theta * (y - lam.ct) / (lam.ct * theta + lam.ct ^ 2)
    )

    d.th = ifelse(
      y == 0,
      denom.y0 * (1 - psi) * (theta / (theta + lam.ct)) ^ theta *
        (log(theta) - log(lam.ct + theta) + lam.ct / (lam.ct + theta)),
      digamma(y + theta) - digamma(theta) + 1 + log(theta) -
        log(theta + lam.ct) - (theta + y) / (theta + lam.ct)
    )

    d.zi = z * d.zi * weights
    d.ct = X * d.ct * weights
    d.th = theta * d.th * weights # multiply by theta due to exponential transform

    return(-unname(colSums(cbind(d.ct, d.zi, d.th))))

  }


  dist      = match_arg(dist, c("pois", "nbinom"))
  link.ct   = match_arg(link.ct, c("sqrt", "identity", "log"))
  link.zi   = match_arg(link.zi, c("logit", "probit", "cauchit", "log", "cloglog"))
  optimizer = match_arg(optimizer, c("nlm", "optim"))


  if (!env_has(e.check, "zi.dens.checks"))
    env_poke(e.check, "zi.dens.checks", T)
  on.exit(rm("zi.dens.checks", envir = e.check), add = T)

  environment(zic.ll) = environment()
  environment(nbinom.grad) = environment()
  environment(pois.grad) = environment()

  nb = dist == "nbinom"

  if (!is.null(fmla)) {
    if (missing(data)) {
      data = environment(fmla)
      if (warn.parent) {
        warning("Data not supplied to function, looking in parent environment.")
      }
    }

    fmla = as.Formula(fmla)

    #get call for model frame
    mf = match.call(expand.dots = F)
    m = match(c("data", "weights", "subset", "na.action"), names(mf), 0)
    mf = mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf$formula = fmla
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    #set up data
    y = model.response(mf)
    X = model.matrix(fmla, mf, rhs = 1)
    z = model.matrix(fmla, mf, rhs = 2)



    weights = if (!is.null(model.weights(mf)))
      model.weights(mf)
    else
      rep(1, length(y))
    offset.ct = get.offset(attr(fmla, "rhs")[[1]], mf)
    offset.zi = get.offset(attr(fmla, "rhs")[[2]], mf)

  }

  assert(
    is.null(starts) ||
      ncol(X) + ncol(z) + sum(grepl("nb", dist)) == length(starts),
    "Vector of starting values (`starts`) is not of correct length."
  )

  #get starting values using glm, or simplex run
  if (is.null(starts)) {
    start.zi = glm.fit(
      z,
      as.integer(y == 0),
      offset = offset.zi,
      weights = weights,
      family = binomial(link = link.zi)
    )$coefficients
    start.ct = glm.fit(
      X[y > 0,],
      y[y > 0],
      offset = offset.ct[y > 0],
      weights = weights[y > 0],
      family = poisson(link = link.ct)
    )$coefficients
    starts = if (dist == "pois")
      c(start.ct, start.zi)
    else
      c(start.ct, start.zi, 0)


  }


  # Check and replace required defaults
  opts = list(...)

  if (!is.null(opts$hessian) && opts$hessian == F)
    warning("Arg `hessian` must be true; changing this automatically.")
  if ((!is.null(opts$control[["fnscale"]]) &&
      opts$control[["fnscale"]] < 0) ||
      (!is.null(opts$fscale) && opts$fscale < 0))
    warning("Function must not be inverted; changing automatically.")
  opts$hessian = T


  # Add in optional defaults
  default.names = switch(
    optimizer,
    "optim" = c("maxit", "reltol"),
    "nlm" = c("iterlim", "gradtol")
  )


  #main optimization, passing additional args via ... using defaults
  res = #
    if (optimizer == "optim") {
      add = list(
        par = starts,
        fn = zic.ll,
        gr = get(paste0(dist, ".grad")),
        grad = F
      )
      opts = c(opts, add)
      opts$control = set.defaults(opts$control, "maxit", 25000)
      opts = set.defaults(opts, "method", "BFGS")
      opts$control[["fnscale"]] = 1
      do.call(optim, opts)
    } else if (optimizer == "nlm") {
      add = list(p = starts,
                 f = zic.ll,
                 grad = T)
      opts = c(opts, add)
      opts = set.defaults(opts, c("iterlim", "gradtol"), c(25000, 1e-8))
      opts$fscale = 1
      alter.cond(do.call(nlm, opts), suppress=T)
    }

  ### clean up results for printing
  coefs = switch(optimizer,
                 "optim" = res$par,
                 "nlm" = res$estimate)
  loglik =  switch(optimizer,
                   "optim" = -res$value,
                   "nlm" = -res$minimum)
  convergence = switch(optimizer,
                       "optim" = res$convergence,
                       "nlm" = res$code)
  niter = switch(optimizer,
                 "optim" = res$counts[1],
                 "nlm" = res$iterations)



  if ((optimizer == "optim" && convergence != 0) ||
      (optimizer == "nlm" && convergence > 2))
    warning(
      paste0(
        "Model did not converge! See `?",
        optimizer,
        "` for more details.\nConvergence code: ",
        convergence
      )
    )

  npar = length(coefs)
  hess = hess.orig = res$hessian
  tol = .Machine$double.eps

  covmat = covmat.orig = tryCatch(
    solve(hess, tol = tol),
    error = function(e)
      matrix(NA, npar, npar)
  )



  #transform hessian if negbin to get correct SEs for dispersion
  if (nb) {
    coefs[npar] = exp(coefs[npar])
    mult = 1 / coefs[npar]

    if (mult < 5e-2) {
      tol = 1e-100
      warning(
        "Dispersion parameter suspiciously large; negative binomial distribution may be inappropriate."
      )
    }

    hess[, npar] = hess[, npar] * mult
    hess[npar, ] = hess[npar, ] * mult
    covmat = solve(hess, tol=tol)

  }


  se = sqrt(diag(covmat))
  ztest = coefs / se
  pval = 2 * pnorm(abs(ztest), lower.tail = F)
  all = list(coefs, se, ztest, pval)
  coefmat.all = do.call(cbind, all)
  aic = 2 * (npar - loglik)
  bic = npar * log(length(y)) - 2 * loglik
  k.ct = ncol(X)
  k.zi = ncol(z)
  if (nb) {
    theta = matrix(sapply(all, tail, 1), nrow = 1)
    rownames(theta) = "Theta"
    colnames(theta) = c("Estimate", "Std. Error", "z-score", "p-value")
    coefmat.ndisp = coefmat.all[-(k.ct + k.zi + 1),]
  } else {
    coefmat.ndisp = coefmat.all
  }

  sum.mat.ct = head(coefmat.ndisp, k.ct)
  sum.mat.zi = tail(coefmat.ndisp, k.zi)
  names(se) = names(coefs) = #
    paste0(c(rep("ct_", k.ct), rep("zi_", k.zi), if (nb)
      ""
      else
        NULL),
      c(colnames(X), colnames(z), if (nb)
        "Theta"
        else
          NULL))

  colnames(sum.mat.ct) = colnames(sum.mat.zi) = colnames(coefmat.all) = c("Estimate",
                                                  "Std. Error",
                                                  "z-score",
                                                  "p-value")
  rownames(sum.mat.ct) = colnames(X)
  rownames(sum.mat.zi) = colnames(z)
  rownames(coefmat.all) = c(colnames(X), colnames(z), if(nb) rownames(theta) else NULL)

  out = list(
    call = match.call(),
    obj = "Zero Inflated Count Model",
    coef = coefs,
    se   = se,
    theta = if (nb)
      theta
    else
      NULL,
    covmat = covmat,
    coefmat.all = coefmat.all,
    coefmat.ct = sum.mat.ct,
    coefmat.zi = sum.mat.zi,
    grad = if (optimizer == "nlm")
      res$gradient
    else
      NULL,
    link.ct = link.ct,
    link.zi = link.zi,
    dist = dist,
    optimizer = optimizer,
    nobs = nrow(X),
    aic = aic,
    bic = bic,
    loglik = loglik,
    convergence = convergence,
    model = if (keep)
      list(
        y = y,
        X = X,
        z = z,
        offset.ct = offset.ct,
        offset.zi = offset.zi,
        weights = weights
      )
    else
      NULL
  )

  class(out) = "zicreg"
  return(out)
}
