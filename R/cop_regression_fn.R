# function to help with getting two offsets instead of one from model frame
get.offset = function(formula, model.frame) {
  name = grep("offset",
              unlist(strsplit(as.character(formula), " \\+ ")),
              value = T)
  
  if (length(name) == 0)
    return(rep(0, nrow(model.frame)))
  else
    return(as.vector(model.frame[, name]))
}

# takes vector of all estimated parameters and subets into coefmat
# to extract each part of each equation (eg. count, zi, eq1, eq2)
make.coefmat = function(beta, se, z, pval, index) {
  if (is.null(index))
    return(NULL)
  else
    out = do.call(cbind, lapply(list(beta, se, z, pval), "[", index))
  
  colnames(out) = c("Estimate", "Std. Err.", "Z value", "Pr(>|z|)")
  return(out)
}


starts.marginal = function() {
  oldwarn = getOption("warn")
  options(warn = -1)
  on.exit(options(warn = oldwarn), add = T)
  
  if (!env_has(e.check, "univ.check")) {
    env_poke(e.check, "univ.check", T)
    on.exit(rm("univ.check", envir = e.check), add = T)
  }
  
  starts.ct = starts.zi = starts.disp = list()
  
  for (i in 1:2) {
    i = as.numeric(i)
    mod <- switch(
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
        zic.reg(y=y[,i], X=X[[i]], z=Z[[i]], 
                offset.ct = offset.ct[[i]], 
                offset.zi=offset.zi[[i]],
          weights = weights,
          dist = "nbinom",
          link.zi = link.zi[i],
          link.ct = link.ct[i],
          gradtol = 1e-4
        ),
      "zip"    =
        zic.reg(
          y=y[,i], X=X[[i]], z=Z[[i]], 
          offset.ct = offset.ct[[i]], 
          offset.zi=offset.zi[[i]],
          dist = "pois",
          link.zi = link.zi[i],
          link.ct = link.ct[i],
          gradtol = 1e-4
        )
    )
    
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
}


# main copula regression function
bizicount = function(fmla1,
                    fmla2,
                    data,
                    cop = "gaus",
                    margins = c("pois", "pois"),
                    link.ct = c("log", "log"),
                    link.zi = c("logit", "logit"),
                    starts = NULL,
                    keep = F,
                    subset,
                    na.action,
                    weights,
                    frech.min = 1e-7,
                    pmf.min = 1e-7,
                    ...) {
  #some arg checking
  check_biv_args()
  
  cop     = match_arg(cop, c("frank", "gaus"))
  link.ct = match_arg(link.ct, c("sqrt", "identity", "log"), several.ok = T)
  link.zi = match_arg(link.zi, c("logit", "probit", "cauchit", "log", "cloglog"), several.ok = T)
  
  # Redefine copula likelihood's environment
  environment(cop.lik) = environment()
  environment(starts.marginal) = environment()
  
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
  kz = unlist(kz)
  if (n.zi == 1)
    zi.ind[[zipos]] = (sum(kx, 1)):(sum(kx, kz))
  if (n.zi == 2) {
    zi.ind[[1]] = (sum(kx, 1)):(sum(kx, kz[1]))
    zi.ind[[2]] = (sum(kx, kz[1], 1)):(sum(kx, kz))
  }

  # Get starting values
  if (is.null(starts))
    starts = starts.marginal()
  
  # Prevent arg checking in CDF to speed up likelihood evaluations by using counter
  # in e.check (environment.check)
  # (checking is done at beginning of this function, so nested checks are redundant)
  if (!env_has(e.check, "zi.dens.checks")){
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
  
  reg.inputs = c(list(f=cop.lik, p = starts), varargs)
  
  ### main regression
  out = withCallingHandlers(
    warning = function(w)
      if (grepl("NA|NaN", w))
        invokeRestart("muffleWarning"),
    do.call(nlm, reg.inputs)
  )
  
  # use numDeriv to get hessian if there are issues with NLM's hessian, provided that the gradient is small, there were iterations
  # and neither the hessian nor gradient are exactly zero (in those situations, we don't want to get a hessian as there are bigger problems)
    if(out$iterations > 1 && 
     all(abs(out$gradient)  < .05) &&
     !all(out$hessian == 0) && 
     !all(out$gradient == 0) &&
     anyNA(suppressWarnings(sqrt(diag(solve(out$hessian)))))
     ){
    
    warning("nlm() was unable to obtain Hessian matrix, so numDeriv::hessian() was used in computing standard errors. 
            Consider reducing 'stepmax' option to nlm to prevent this.
            See `?nlm` for more details on the 'stepmax' option.")
    out$hessian = numDeriv::hessian(cop.lik, out$estimate, method.args = list(r=6))
  
  }
     
  conv.message = "try adjusting stepmax. See '?nlm' Details --> Value --> Code for more information."
  switch(out$code,
           invisible(NULL),
           warning(paste("Convergence code 2", conv.message)),
           warning(paste("Convergence code 3", conv.message)),
           stop("Maximum iterations reached, increase `iterlim`. IE, add 'iterlim = [some large integer]' to function call."),
           stop(paste("Convergence code 5", conv.message))
  )
  
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
      
      hess.mult[d,] = thpr.inv # replace hess multiplier row d with 1/derivative of dependence
      hess.mult[, d] = hess.mult[, d] * thpr.inv # multiply hess multiplier col d by 1/deriv dependence
      hess.mult[(d - 1),] = hess.mult[(d - 1),] * dpr.inv[1] # multiply second to last row of hess by 1/deriv disp.2 (recall use of rev())
      hess.mult[, (d - 1)] = hess.mult[, (d - 1)] * dpr.inv[1]
      
      if (n.nb > 1) {
        hess.mult[(d - 2) ,] = hess.mult[(d - 2) ,] * dpr.inv[2] #first dispersion transformation (recall use of rev())
        hess.mult[, (d - 2)] = hess.mult[, (d - 2)] * dpr.inv[2]
      }
    }
    
    #transform original hess element-wise
    hess.new = out$hessian * hess.mult
  } else
    hess.new = out$hessian
  
  
  # check that hessian of loglik is neg def
  evs = eigen(-hess.new, symmetric = T)$values
  assert(all(evs <= sqrt(.Machine$double.eps) * abs(evs[1])),
         "Hessian of loglik is not negative definite at convergence point; convergence point is not a maximum.",
         type="warning"
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
      margins = lapply(margins, function(x) x)
    ),
    mapply(
      function(names, margins)
        if (grepl("nb", margins))
          paste0("disp_", names)
      else
        NULL,
      names = colnames(y),
      margins = lapply(margins, function(x) x)
    ),
    "dependence"
  ))
  
  
  # more informative name must be made AFTER coefmats
  # due to print method for output object
  append = c(rep("ct1_", kx[1]),
             rep("ct2_", kx[2]),
             if (n.zi == 1)
               rep(paste0("zi", grep("zi", margins), "_"), kz),
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
        NULL
    )
  
  
  class(res) = "bizicount"
  return(res)
  
}


