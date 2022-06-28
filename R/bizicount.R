# main copula regression function

#' @name bizicount
#' @title Bizicount: Maximum likelihood estimation of copula-based bivariate zero-inflated
#'   (and non-inflated) count models
#'
#' @description The main bivariate regression function of the \code{\link{bizicount-package}}
#'   Estimates copula-based bivariate zero-inflated (and non-inflated)
#'   count models via maximum likelihood. Supports the Frank and Gaussian
#'   copulas, as well as zero-inflated Poisson and negative binomial margins
#'   (and their non-inflated counterparts). It's class has associated
#'   \code{\link[stats]{simulate}} methods for post-estimation diagnostics using
#'   the \code{\link[=DHARMa]{DHARMa}} package, as well as an
#'   \code{\link[texreg]{extract}} method for printing professional tables using
#'   \code{\link[texreg]{texreg}}. See the 'See Also' section for links to these
#'   methods.
#'
#' @details
#' \itemize{
#'  \item \code{starts} -- Starting values should be organized as
#'  follows:
#'    \enumerate{
#'      \item count parameters for margin 1
#'      \item count parameters for margin 2
#'      \item zero-inflated parameters for margin 1 (if applicable),
#'      \item zero-inflated parameters for margin 2 (if applicable),
#'      \item inverse dispersion parameter for margin 1 (if applicable),
#'      \item inverse dispersion parameter for margin 2 (if applicable)
#'    }
#'  Thus, in general count parameters should come first, followed by
#'  zero-inflation parameters, and finally inverse dispersion parameters.
#'
#' \item \code{frech.min} -- Changing this argument should almost never be
#' necessary. Frechet (1951) and Hoeffding (1940) showed that copula CDFs have
#' bounds of the form \eqn{max{u + v - 1, 0} \le C(u, v) \le min{u, v}}, where
#' \eqn{u} and \eqn{v} are uniform realizations derived from the probability
#' integral transform. Due to numerical underflow, very small values of \eqn{u}
#' and \eqn{v} can be rounded to zero. Particularly when evaluating the Gaussian
#' copula CDF this is problematic, ultimately leading to infinite-valued
#' likelihood evaluations. Therefore, we impose Frechet-Hoeffding bounds
#' numerically as \eqn{max{u + v - 1, frech.min} \le C(u, v) \le min{u, v, 1 -
#' frech.min}}. NOTE: Setting this to 0 imposes the original Frechet bounds
#' mentioned above.
#'
#' \item \code{pmf.min} -- Changing this argument should almost never be
#' necessary. Observations can have likelihoods that are extremely close to 0.
#' Numerically, these get rounded to 0 due to underflow. Then, taking logarithms
#' results in an infinite likelihood. To avoid this, we bound PMF evaluations
#' from below at \code{pmf.min}.
#'
#' \item `...` -- Sometimes it may be useful to alter \code{\link[stats]{nlm}}'s
#' default parameters. This can be done by simply passing those arguments into
#' `bizicount()`. The two that seem to benefit the fitting process the most are
#' `stepmax` and `steptol`. Readers are referred to the documentation on
#' \code{\link[stats]{nlm}} for more details on these parameters. It can be
#' useful to lower `stepmax` particularly when the Hessian is not negative
#' definite at convergence, sometimes to a value between 0 and 1. It can also be
#' beneficial to increase `steptol`.
#'
#' }
#'
#'
#'
#' @references Genest C, Nešlehová J (2007). “A primer on copulas for count
#'   data.” ASTIN Bulletin: The Journal of the IAA, 37(2), 475–515.
#'
#'   Inouye DI, Yang E, Allen GI, Ravikumar P (2017). “A review of multivariate
#'   distributions for count data derived from the Poisson distribution.” Wiley
#'   Interdisciplinary Reviews: Computational Statistics, 9(3).
#'
#'   Joe H (1997). Multivariate models and multivariate dependence concepts. CRC Press.
#'
#'   Nikoloulopoulos A (2013). “Copula-Based Models for Multivariate Discrete
#'   Response Data.” In P Jaworski, F Durante, WK Härdle (eds.), Copulae in
#'   Mathematical and Quantitative Finance, chapter 11, pp. 231–250. Springer.
#'
#'   Nelsen RB (2007). An Introduction to Copulas. Springer Science & Business Media.
#'
#'   Trivedi P, Zimmer D (2017). “A note on identification of bivariate copulas
#'   for discrete countdata.” Econometrics, 5(1), 10.
#'
#'   Trivedi PK, Zimmer DM (2007). Copula modeling: an introduction for
#'   practitioners. NowPublishers Inc.
#'
#' @example inst/examples/bizicount_ex.R
#'
#' @return An S3 \code{\link{bizicount-class}} object, which is a list containing:
#' \itemize{
#' \item `coef` -- Coefficients of the model
#' \item `coef.nid` -- Coefficients without margin IDs
#' \item `coef.orig` -- Coefficients prior to transformations, for Gaussian
#'   dependence and negative binomial dispersion.
#' \item `coef.orig.nid` -- Coefficients prior to transforms, no margin IDs.
#' \item `se` -- Asymptotic normal-theory standard errors based on observed Fisher Information
#' \item `se.nid` -- Standard errors without margin IDs
#' \item `z` -- z-scores for parameter estimates
#' \item `z.nid` -- z-scores without margin IDs
#' \item `p` -- p-values for parameter estimates
#' \item `p.nid` -- p-values without margin IDs
#' \item `coefmats` -- A list containing coeficient matrices for each margin
#' \item `loglik` -- Scalar log-likelihood at convergence
#' \item `grad` -- Numerical gradient vector at convergence
#' \item `n.iter` -- Number of quasi-newton fitting iterations.
#' \item `covmat` -- Covariance matrix of parameter estimates based on observed Fisher Information
#' \item `aic` -- Model's Akaike information
#' \item `bic` -- Model's Bayesian information criterion
#' \item `nobs` -- Number of observations
#' \item `margins` -- Marginal distributions used in fitting
#' \item `link.zi, link.ct` -- Names of link functions used in fitting
#' \item `invlink.ct, invlink.zi` -- Inverse link functions used in fitting (the
#'   actual function, not their names)
#' \item `outcomes` -- Name of the response vector
#' \item `conv` -- Integer telling convergence status in nlm. See ?nlm.
#' \item `cop` -- The copula used in fitting
#' \item `starts` -- list of starting values used
#' \item `call` -- The model's call
#' \item `model` -- List containing model matrices, or `NULL` if `keep = F`.
#' \item `scaled` -- Vector indicating which covariates in each margin were scaled
#' according to the `scaling` parameter.
#' }
#'
#'
#'
#' @param fmla1,fmla2 \code{\link[stats]{formula}}s for the first margin and
#'   second margins, respectively. If non-inflated, of the form `y ~ x_1 + x_2 +
#'   ... + x_k`; if inflated, of the form `y ~ x1 + x2 + ... + x_k| z1 + z2 +
#'   ... + z_p`, where `y` is the outcome for the first margin, `x` are
#'   covariates for count parameters, and `z` are covariates for zero-inflated
#'   parameters in each margin. All covariates can be the same.
#' @param data A \code{\link[base]{data.frame}} containing the response variables, covariates, and
#'   offsets for the model. If `NULL`, these quantities are searched for in the
#'   parent environment.
#' @param cop Character string specifying the copula to be used. One of
#'   `c("gaus", "frank")`. Partial matching supported.
#' @param margins Length 2 character vector specifying the marginal
#'   distributions for each outcome. Each of the two elements must be one of
#'   `c("pois", "nbinom", "zip", "zinb")`, and must be consistent with its
#'   corresponding formula (i.e., zero-inflated margins with zero-inflated
#'   formulas).
#' @param link.ct Length 2 character string specifying the link function used
#'   for the count portion of each margin. One of `c("log", "identity",
#'   "sqrt")`.
#' @param link.zi Length 2 character string specifying the link function used
#'   for the zero-inflation portion of each margin. One of `c("logit", "probit",
#'   "cauchit", "log", "cloglog")`. Ignored if corresponding `margins` entry is
#'   not zero-inflated.
#' @param scaling Character string (partial matching supported) indicating what
#' type of scaling to apply to covariates (if any). One of `c("none", "1sd", "gelman", "mm")`.
#'    * `"none"` will apply no alterations to covariates.
#'    * `"1sd"` will subtract off the mean, and divide by one standard deviation.
#'    * `"gelman"` will subtract off the mean, and divide by two standard deviations,
#'       which makes binary and continuous variables have a similar interpretation.
#'       See Gelman (2008), "Scaling Regression Inputs by Dividing by Two Standard Deviations."
#'    * `"mm"` will apply min-max normalization so that continuous covariates lie within a unit hypercube.
#'
#'  Factor variables, offsets, and the intercept are not scaled. If scaling,
#'  it is recommended that data be supplied in a dataframe (as opposed to from the global environment),
#'  otherwise the automated scaling cannot reliably be applied. The names of variables that have been
#'  scaled are returned as part of the `bizicount` object, in the list-element called `scaled`.
#' @param starts Numeric vector of starting values for parameter estimates. See
#'   'Details' section regarding the correct order for the values in this vector.
#'   If `NULL`, starting values are obtained automatically by a univariate regression fit.
#' @param keep Logical indicating whether to keep the model matrix in the
#'   returned model object. Defaults to `FALSE` to conserve memory. NOTE: This
#'   must be set to `TRUE` to use the  \code{\link[texreg]{texreg}},
#'   \code{\link{simulate.bizicount}}, \code{\link[stats]{fitted}}, or
#'   \code{\link{make_DHARMa}} functions with `bizicount` objects.
#' @param subset A vector indicating the subset of observations to use in
#' estimation.
#' @param na.action A function which indicates what should happen when the data
#'   contain NAs. Default is \code{\link[stats]{na.omit}}.
#' @param weights An optional numeric vector of weights for each observation.
#' @param frech.min Lower boundary for Frechet-Hoeffding bounds on copula CDF.
#'   Used for computational purposes to prevent over/underflow in likelihood
#'   search. Must be in \eqn{[0, 1e-5]}, with \eqn{0} imposing the original FH
#'   bounds without computational consideration. See 'Details.'
#' @param pmf.min Lower boundary on copula PMF evaluations. Used for
#'   computational purposes to prevent over/underflow in likelihood search. Must
#'   be in \eqn{[0, 1e-5]}, with \eqn{0} imposing no bound. See `Details.'
#' @param ... Additional arguments to be passed on to the quasi-newton fitting
#'   function, \code{\link[stats]{nlm}}. See 'Details' for some parameters that
#'   may be useful to alter.
#'
#' @author John Niehaus
#'
#'
#' @export
bizicount = function(fmla1,
                     fmla2,
                     data,
                     cop = "gaus",
                     margins = c("pois", "pois"),
                     link.ct = c("log", "log"),
                     link.zi = c("logit", "logit"),
                     scaling = "none",
                     starts = NULL,
                     keep = FALSE,
                     subset,
                     na.action,
                     weights,
                     frech.min = 1e-7,
                     pmf.min = 1e-7,
                     ...) {
  #some arg checking
  check_biv_args()

  # Bivariate likelihood (univariate is lower down)
  cop.lik = function(parms) {
    ct = lapply(ct.ind, function(x)
      parms[x])


    # Loop over inputs to create list of marginal lambdas
    # ie, lam[1] = f1(x1, z1, w1), lam[2] = f2(x2, z2, w2)
    lam = mapply(
      function(par, design, offset, invlink) {
        if (ncol(design) == 1)
          as.vector(invlink(design * par + offset))
        else
          as.vector(invlink(design %*% par + offset))
      },
      par = ct,
      design = X,
      offset = offset.ct,
      invlink = invlink.ct,
      SIMPLIFY = F
    )

    # Get zi parameters via subset
    zi = list()
    if (n.zi == 1)
      zi[[zipos]] = parms[zi.ind[[zipos]]]
    if (n.zi == 2) {
      zi = lapply(zi.ind, function(x)
        parms[x])
    }

    # dispersion and dependence params
    disp = list()
    if (n.nb == 1)
      disp[[nbpos]] = exp(tail(parms, 2)[1])
    if (n.nb == 2) {
      disp[[1]] = exp(tail(parms, 3)[1])
      disp[[2]] = exp(tail(parms, 2)[1])
    }

    dep = tail(parms, 1)
    if (cop == 'gaus')
      dep = tanh(dep)

    # Create list of zero-inflation probability vectors
    psi = list()
    if (n.zi > 0) {
      for (i in zipos)
        psi[[i]] = as.vector(invlink.zi[[i]](Z[[i]] %*% zi[[i]] + offset.zi[[i]]))
    }


    # get list of arguments to marginal cdfs
    margin.args = list()

    for (i in c(1, 2)) {
      margin.args[[i]] = switch(
        margins[i],
        "pois" = list(q      = y[, i],
                      lambda = lam[[i]]),
        "nbinom"  = list(
          q     = y[, i],
          mu    = lam[[i]],
          size  = disp[[i]]
        ),
        "zip" = list(
          q      = y[, i],
          lambda = lam[[i]],
          psi    = psi[[i]]
        ),
        "zinb" = list(
          q     = y[, i],
          mu    = lam[[i]],
          size  = disp[[i]],
          psi   = psi[[i]]
        )
      )
    }



    # evaluate copula pmf
    d = cop.pmf(
      margin.args = margin.args,
      margins = margins,
      cop = cop,
      dep = dep,
      pmf.min = pmf.min,
      frech.min = frech.min
    )

    # bound likelihoods for logging
    d = bound(d * weights, low = pmf.min)


    #ll

    -sum(log(d))

  }

  starts.marginal = function() {
    if (!env_has(e.check, "univ.check")) {
      env_poke(e.check, "univ.check", T)
      on.exit(rm("univ.check", envir = e.check), add = T)
    }

    starts.ct = starts.zi = starts.disp = list()

    for (i in 1:2) {
      i = as.numeric(i)
      mod <- suppressWarnings(switch(
        margins[i],
        "pois"   =
          glm.fit(
            X[[i]],
            y[, i],
            offset = offset.ct[[i]],
            family = poisson(link = link.ct[i]),
            control = list(epsilon = 1e-4),
            weights = weights
          ),
        "nbinom" = suppressWarnings(eval.parent(substitute(
          MASS::glm.nb(
            y[, i] ~ X[[i]] + 0 + offset(offset.ct[[i]]),
            weights = weights,
            link = link.ct[i]
          ),
          list(
            link.ct = get("link.ct", parent.frame()),
            i = i
          )
        ))),
        "zinb"   =
          zic.reg(
            y = y[, i],
            X = X[[i]],
            z = Z[[i]],
            offset.ct = offset.ct[[i]],
            offset.zi = offset.zi[[i]],
            weights = weights,
            dist = "nbinom",
            link.zi = link.zi[i],
            link.ct = link.ct[i],
            gradtol = 1e-4
          ),
        "zip"    =
          zic.reg(
            y = y[, i],
            X = X[[i]],
            z = Z[[i]],
            offset.ct = offset.ct[[i]],
            offset.zi = offset.zi[[i]],
            dist = "pois",
            link.zi = link.zi[i],
            link.ct = link.ct[i],
            gradtol = 1e-4
          )
      ))

      if (grepl("nb", margins[i]))
        starts.disp[[i]] = mod$theta[1]

      coefs = coef(mod)

      starts.ct[[i]] = coefs[1:kx[i]]
      starts.zi[[i]] = if (any(zipos == i))
        coefs[(kx[i] + 1):(kx[i] + kz[[i]])]
      else
        NULL
    }

    return(na.omit(c(
      unlist(starts.ct),
      unlist(starts.zi),
      unlist(starts.disp),
      cor(y, method = "spearman")[1, 2]
    )))
  } # end of starting values function

  scaling = match_arg(scaling, c("none", "1sd", "gelman", "mm"))
  scaling = switch(scaling, "none" = 0, "1sd" = 1, "gelman" = 2, "mm" = 3)
  cop     = match_arg(cop, c("frank", "gaus"))
  link.ct = match_arg(link.ct, c("sqrt", "identity", "log"), several.ok = T)
  link.zi = match_arg(link.zi,
                      c("logit", "probit", "cauchit", "log", "cloglog"),
                      several.ok = T)


  # indices for getting parms/matrices etc from formulas
  n.zi = sum(grepl("zi", margins))
  l.zi = any(grepl("zi", margins))
  zipos = grep("zi", margins)

  n.nb = sum(grepl("nb", margins))
  l.nb = any(grepl("nb", margins))
  nbpos = grep("nb", margins)

  if (missing(data)) {
    warning("Data not supplied to function, looking in parent environment.",
            immediate. = T)
    data = environment(fmla1)
  }

  #get model matrices, offsets, weights
  fmla = as.Formula(fmla1, fmla2)
  fmla.list = list(as.Formula(fmla1),
                   as.Formula(fmla2))

  mf = match.call(expand.dots = F)
  m = match(c("data", "weights", "subset", "na.action"), names(mf), 0)
  mf = mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula = fmla
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  y = model.part(fmla, mf, lhs = c(1, 2))

  X = lapply(fmla.list, function(x)
    model.matrix(x, mf, rhs = 1))


  Z = list()
  for (i in zipos)
    Z[[i]] = model.matrix(fmla.list[[i]], mf, rhs = 2)

  if (length(Z) == 1)
    Z[[2]] = NULL


  if (scaling != 0){
    X = lapply(X, scaler, scaling = scaling)
    Z = lapply(Z, scaler, scaling = scaling)
  }

  offset.ct = lapply(fmla.list,
                     function(x)
                       get.offset(attr(x, "rhs")[[1]], mf))

  offset.zi = list()
  for (i in zipos)
    offset.zi[[i]] = get.offset(attr(fmla.list[[i]], "rhs")[[2]], mf)

  weights = #
    if (!is.null(model.weights(mf)))
      model.weights(mf)
  else
    rep(1, nrow(y))

  #--end model matrix stuff

  # List of link functions, num of params in each part
  invlink.ct = lapply(link.ct, function(x)
    make.link(x)$linkinv)
  invlink.zi = lapply(link.zi, function(x)
    make.link(x)$linkinv)
  kx = sapply(X, ncol)
  kz = lapply(Z, ncol)
  # Get parameter indices corresponding to each part of each margin
  ct.ind = list(1:kx[1],
                (kx[1] + 1):(sum(kx)))


  zi.ind = list()
  if (n.zi == 1)
    zi.ind[[zipos]] = (sum(kx, 1)):(sum(kx, unlist(kz)))
  if (n.zi == 2) {
    zi.ind[[1]] = (sum(kx, 1)):(sum(kx, kz[[1]]))
    zi.ind[[2]] = (sum(kx, kz[[1]], 1)):(sum(kx, unlist(kz)))
  }

  # Get starting values
  if (is.null(starts))
    starts = starts.marginal()


  # Prevent arg checking in CDF to speed up likelihood evaluations by using counter
  # in e.check (environment.check)
  # (checking is done at beginning of this function, so nested checks are redundant)
  if (!env_has(e.check, "zi.dens.checks")) {
    env_poke(e.check, "zi.dens.checks", T)
    on.exit(rm("zi.dens.checks", envir = e.check), add = T)
  }

  # Additional args to pass on to optimization
  varargs = list(...)

  # Set defaults
  if (is.null(varargs$iterlim))
    varargs$iterlim = 100000
  if (!is.null(varargs$hessian) && !isTRUE(varargs$hessian))
    warning("Arg `hessian` must be TRUE. Changing automatically...",
            immediate. = T)

  varargs$hessian = T

  reg.inputs = c(list(f = cop.lik, p = starts), varargs)

  ### main regression
  out = withCallingHandlers(
    warning = function(w)
      if (grepl("NA|NaN", w))
        invokeRestart("muffleWarning"),
    do.call(nlm, reg.inputs)
  )

  # use numDeriv to get hessian if there are issues with NLM's hessian, provided that the gradient is small, there were iterations
  # and neither the hessian nor gradient are exactly zero (in those situations, we don't want to get a hessian as there are bigger problems)
  if (out$iterations > 1 &&
      all(abs(out$gradient)  < .05) &&
      !all(out$hessian == 0) &&
      !all(out$gradient == 0) &&
      anyNA(tryCatch(
        suppressWarnings(sqrt(diag(solve(
          out$hessian
        )))),
        error = function(e)
          return(NA)
      ))) {
    warning(
      "nlm() was unable to obtain Hessian matrix, so numDeriv::hessian() was used in computing standard errors.
            Consider reducing 'stepmax' option to nlm to prevent this.
            See `?nlm` for more details on the 'stepmax' option."
    )

    out$hessian = numDeriv::hessian(cop.lik, out$estimate, method.args = list(r = 6))

  }

  conv.message = "try adjusting stepmax. See '?nlm' Details --> Value --> Code for more information."
  switch(out$code,
         invisible(NULL),
         warning(paste("Convergence code 2", conv.message)),
         warning(paste("Convergence code 3", conv.message)),
         stop(
           "Maximum iterations reached, increase `iterlim`. IE, add 'iterlim = [some large integer]' to function call."
         ),
         stop(paste("Convergence code 5", conv.message)))

  # transform hessian with change of variables to get correct standard errors for
  # transformed parameters. only necessary if negative binomial margins or gaussian copula
  if (n.nb > 0 || cop == "gaus") {
    hess.mult = matrix(1,
                       nrow = nrow(out$hessian),
                       ncol = ncol(out$hessian))


    thpr.inv = if (cop == "gaus")
      cosh(tail(out$estimate, 1)) ^ 2
    else
      1

    if (n.nb > 0) {
      # get dispersion, 1/derivative, reverse them for multiplication below
      dpr.inv = rev(1 / exp(tail(out$estimate, (n.nb + 1))[1:n.nb]))
      d = nrow(hess.mult)

      hess.mult[d, ] = thpr.inv # replace hess multiplier row d with 1/derivative of dependence
      hess.mult[, d] = hess.mult[, d] * thpr.inv # multiply hess multiplier col d by 1/deriv dependence
      hess.mult[(d - 1), ] = hess.mult[(d - 1), ] * dpr.inv[1] # multiply second to last row of hess by 1/deriv disp.2 (recall use of rev())
      hess.mult[, (d - 1)] = hess.mult[, (d - 1)] * dpr.inv[1]

      if (n.nb > 1) {
        hess.mult[(d - 2) , ] = hess.mult[(d - 2) , ] * dpr.inv[2] #first dispersion transformation (recall use of rev())
        hess.mult[, (d - 2)] = hess.mult[, (d - 2)] * dpr.inv[2]
      }
    }

    #transform original hess element-wise
    hess.new = out$hessian * hess.mult
  } else
    hess.new = out$hessian


  # check that hessian of loglik is neg def
  evs = eigen(-hess.new, symmetric = T)$values
  assert(
    all(evs <= sqrt(.Machine$double.eps) * abs(evs[1])),
    "Hessian of loglik is not negative definite at convergence point; convergence point is not a maximum.",
    type = "warning"
  )

  #cleanup results for print
  beta = beta.orig = out$estimate

  npar = length(beta)

  if (n.nb == 1)
    beta[(npar - 1)] = exp(beta[(npar - 1)])
  if (n.nb == 2)
    beta[(npar - 2):(npar - 1)] = exp(tail(beta, 3)[1:2])
  if (cop == "gaus")
    beta[npar] = tanh(tail(beta, 1))

  if (n.nb > 0 && any(dpr.inv < 5e-2)) {
    tol = 1e-100
    warning(
      "Dispersion parameter is suspiciously large; negative binomial margins may be inappropriate."
    )
  } else {
    tol = .Machine$double.eps
  }

  se = tryCatch(
    sqrt(diag(solve(hess.new, tol = tol))),
    error = function(e)
      rep(NA, npar)
  )

  z  = beta / se
  p  = 2 * pnorm(abs(z), lower.tail = F)

  dep.ind = length(beta)
  disp.ind = if (n.nb > 0)
    tail((1:npar), (n.nb + 1))[1:n.nb]
  else
    NULL

  #names
  namevec = c(
    unlist(sapply(X, colnames)),
    unlist(sapply(Z, colnames)),
    if (n.nb > 0)
      paste0("disp.", grep("nb", margins))
    else
      NULL,
    "dependence"
  )

  names(beta) = names(beta.orig) = names(se) = names(z) = names(p) = namevec
  colnames(hess.new) = colnames(out$hessian) = #
    rownames(hess.new) = rownames(out$hessian) = namevec

  #get coef matrices
  coefmats = list()
  all.inds = c(ct.ind, zi.ind, disp.ind, dep.ind)
  coefmats = lapply(
    all.inds,
    make.coefmat,
    beta = beta,
    se = se,
    z = z,
    pval = p
  )

  coefmats = coefmats[!sapply(coefmats, is.null)]

  names(coefmats) = unlist(c(
    paste0("ct_", colnames(y)),
    mapply(
      function(names, margins)
        if (grepl("zi", margins))
          paste0("zi_", names)
      else
        NULL,
      names = colnames(y),
      margins = lapply(margins, function(x)
        x)
    ),
    mapply(
      function(names, margins)
        if (grepl("nb", margins))
          paste0("disp_", names)
      else
        NULL,
      names = colnames(y),
      margins = lapply(margins, function(x)
        x)
    ),
    "dependence"
  ))


  # more informative name must be made AFTER coefmats
  # due to print method for output object
  append = c(rep("ct1_", kx[1]),
             rep("ct2_", kx[2]),
             if (n.zi == 1)
               rep(paste0("zi", grep("zi", margins), "_"), unlist(kz)),
             if (n.zi == 2) {
               c(rep("zi1_", kz[1]),
                 rep("zi2_", kz[2]))
             })

  # Save coefs, se, etc without equation ids
  beta.nid = beta
  beta.orig.nid = beta.orig
  se.nid = se
  p.nid = p
  z.nid = z

  # Add equation ids to coefs, se, etc
  appind = seq_along(append)
  names(beta)[appind] = names(beta.orig)[appind] = #
    names(se)[appind] = names(z)[appind] = #
    names(p)[appind] = #
    paste0(append, names(beta)[appind])

  # indicate which vars are scaled
  scaled = unique(unlist(c(
    lapply(X, attr, which = "scaled"),
    lapply(Z, attr, which = 'scaled')
    )))


  res =
    list(
      coef = beta,
      coef.nid = beta.nid,
      coef.orig = beta.orig,
      coef.orig.nid = beta.orig.nid,
      se = se,
      se.nid = se.nid,
      z = z,
      z.nid = z.nid,
      p = p,
      p.nid = p.nid,
      coefmats = coefmats,
      loglik = -out$minimum,
      grad = -out$gradient,
      n.iter = out$iterations,
      covmat = suppressWarnings(solve(hess.new, tol = tol)),
      aic = 2 * (npar - (-out$minimum)),
      bic = npar * log(length(y)) - 2 * (-out$minimum),
      nobs = nrow(y),
      margins = margins,
      link.zi = link.zi,
      link.ct = link.ct,
      invlink.ct = invlink.ct,
      invlink.zi = invlink.zi,
      outcomes = colnames(y),
      conv = out$code,
      cop = cop,
      starts = starts,
      call = match.call(expand.dots = T),
      model = if (keep)
        list(
          y = y,
          X = X,
          Z = Z,
          offset.ct = offset.ct,
          offset.zi = offset.zi,
          offset = offset.zi,
          weights = weights
        )
      else
        NULL,
      scaled = scaled
    )


  class(res) = "bizicount"
  return(res)

}
