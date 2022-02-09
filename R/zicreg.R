#univariate regression function
#' @export
#' @name zic.reg
#' @title Univariate zero-inflated Poisson and negative binomial regression models
#'
#' @description This function from the \code{\link{bizicount}} package estimates
#'   univariate zero-inflated Poisson and negative binomial regression models
#'   via maximum likelihood using either the \code{\link[stats]{nlm}} or
#'   \code{\link[stats]{optim}} optimization functions.  It's class has
#'   associated \code{\link[stats]{simulate}} methods for post-estimation
#'   diagnostics using the \code{\link[=DHARMa]{DHARMa}} package, as well as an
#'   \code{\link[texreg]{extract}} method for printing professional tables using
#'   \code{\link[texreg]{texreg}}. Visit the 'See Also' section for links to these
#'   methods for `zicreg` objects.
#'
#' @param fmla A \code{\link[stats]{formula}} of the form `y ~ x_1 + x_2 + ... +
#'   x_n + offset(count_var) | z_1  + ... z_n + offset(zi_var)`, where the `x`
#'   values are covariates in the count portion of the model, and `z` are in the
#'   zero-inflation portion. The `z` and `x` variables can be the same. If `NULL`,
#'   design matrices, the response vector, and offsets can be entered directly; see
#'   `X`, `z`, `y`, `offset.ct`, and `offset.zi` below.
#' @param data A \code{\link[base]{data.frame}} containing all variables
#'   appearing in `fmla`, including offsets. If not specified, variables are
#'   searched for in parent environment.
#' @param dist The distribution used for the count portion of the zero-inflated
#'   mixture. One of `c("pois", "nbinom")`, partial matching supported.
#' @param link.ct String specifying the link function used for the count portion
#'   of the mixture distribution. One of `c("log", "identity", "sqrt")`.
#'   See \code{\link[stats]{family}}.
#' @param link.zi Character string specifying the link function used for the
#'   zero-inflation portion of the mixture distribution. One of `c("logit",
#'   "probit", "cauchit", "log", "cloglog")`. See \code{\link[stats]{family}}.
#' @param optimizer String specifying the optimizer to be used for fitting, one
#'   of `c("nlm", "optim")`. If `"optim"`, defaults to `method="BFGS"`.
#' @param starts Optional vector of starting values used for the numerical
#'   optimization procedure. Should have count parameters first (with intercept
#'   first, if applicable), followed by zero-inflated parameters (with intercept
#'   first, if applicable), and the inverse dispersion parameter last (if
#'   applicable).
#' @param subset Vector indicating the subset of observations on which to
#' estimate the model
#' @param na.action A function which indicates what should happen when the data
#'   contain NAs. Default is \code{\link[stats]{na.omit}}.
#' @param weights An optional numeric vector of weights for each observation.
#' @param X,z If `fmla = NULL`, these are the design matrices of covariates for
#'   the count and zero-inflation portions, respectively. Both require no
#'   missingness. Similar in spirit to \code{\link[stats]{glm.fit}} in that it
#'   can be faster for larger datasets because it bypasses model matrix
#'   creation.
#' @param y If `fmla = NULL`, a vector containing the response variable.
#' @param offset.ct,offset.zi If `fmla = NULL`, vectors containing the
#'   (constant) offset for the count and zero-inflated portions, respectively.
#'   Must be equal in length to `y`, and row-dim of `X`, `z`. If left `NULL`,
#'   defaults to `rep(0, length(y))`.
#' @param warn.parent Logical indicating whether to warn about `data` not
#'   being supplied.
#' @param keep Logical indicating whether to keep the model matrices in the
#'   returned model object. Must be `TRUE` to use \code{\link[=DHARMa]{DHARMa}}
#'   and \code{\link[texreg]{texreg}} with the model object, e.g., via
#'   \code{\link{simulate.zicreg}} and \code{\link{extract.zicreg}}, as well as
#'   base generics like \code{\link[stats]{fitted}} and
#'   \code{\link[stats]{predict}}.
#' @param ... Additional arguments to pass on to the chosen optimizer, either
#' \code{\link[stats]{nlm}} or \code{\link[stats]{optim}}. See 'Examples'.
#'
#' @example /inst/examples/zicreg_ex.R
#'
#' @return An S3 \code{\link{zicreg-class}} object, which is a list containing:
#' \itemize{
#' \item  `call` -- The original function call
#' \item  `obj` -- The class of the object
#' \item  `coef` -- Vector of coefficients, with count, then zi, then dispersion.
#' \item  `se` -- Vector of asymptotic standard errors
#' \item  `grad` -- Gradient vector at convergence
#' \item  `link.ct` -- Name of link used for count portion
#' \item  `link.zi` -- Name of link used for zero-inflated portion
#' \item  `dist` -- Name of distribution used for count portion
#' \item  `optimizer` -- Name of optimization package used in fitting
#' \item  `coefmat.ct` -- Coefficient matrix for count portion
#' \item  `coefmat.zi` -- Coefficient matrix for zero-inflated portion
#' \item  `convergence` -- Convergence code from optimization routine.
#' \item  `coefmat.all` -- Coefficient matrix for both parts of the model
#' \item  `theta` -- Coefficient matrix for dispersion, if applicable.
#' \item  `covmat` -- Asymptotic covariance matrix
#' \item  `nobs` -- Number of observations
#' \item  `aic` -- Akaike information
#' \item  `bic` -- Bayes information
#' \item  `loglik` -- Log-likelihood at convergence
#' \item  `model` -- List containing model matrices if `keep = TRUE`
#'}
#'
#' @author John Niehaus
#'
#' @references Lambert, Diane. "Zero-inflated Poisson regression, with an
#'   application to defects in manufacturing." Technometrics 34.1 (1992): 1-14.
#' @seealso \code{\link{simulate.zicreg}}, \code{\link{extract.zicreg}}

zic.reg = function(fmla = NULL,
                   data,
                   dist = "pois",
                   link.ct = "log",
                   link.zi = "logit",
                   optimizer = "nlm",
                   starts = NULL,
                   subset,
                   na.action,
                   weights = rep(1, length(y)),
                   X = NULL,
                   z = NULL,
                   y = NULL,
                   offset.ct = NULL,
                   offset.zi = NULL,
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
    d.th = theta * d.th * weights # multiply by theta due to exponential transform/chain rule

    return(-unname(colSums(cbind(d.ct, d.zi, d.th))))

  }


  dist      = match_arg(dist, c("pois", "nbinom"))
  link.ct   = match_arg(link.ct, c("sqrt", "identity", "log"))
  link.zi   = match_arg(link.zi, c("logit", "probit", "cauchit", "log", "cloglog"))
  optimizer = match_arg(optimizer, c("nlm", "optim"))


  if (!env_has(e.check, "zi.dens.checks"))
    env_poke(e.check, "zi.dens.checks", T)
  on.exit(rm("zi.dens.checks", envir = e.check), add = T)

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

  if(is.null(offset.ct) && is.null(fmla))
    offset.ct = rep(0, length(y))

  if(is.null(offset.zi) && is.null(fmla))
    offset.zi = rep(0, length(y))

  assert(
    is.null(starts) ||
      ncol(X) + ncol(z) + sum(grepl("nb", dist)) == length(starts),
    "Vector of starting values (`starts`) is not of correct length."
  )

  #get starting values using glm
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
