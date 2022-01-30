#' @importFrom texreg extract
#' @importFrom texreg createTexreg
#' @importFrom DHARMa createDHARMa



#' @export
print.bizicount= function(x, ...){

  for(i in 1:2){
    cat(paste("\nCount:", x$outcomes[i]), "\n")
    print(x$coef[grepl(paste0("ct",i), names(x$coef))])
    if(grepl("zi", x$margins[i])){
      cat(paste("\nZero Infl:", x$outcomes[i]), "\n")
      print(x$coef[grepl(paste0("zi",i), names(x$coef))])
  }}

  disps = sum(grepl("disp", names(x$coef)))

  cat("\n")
  print(tail(x$coef, 1 + disps))

  cat("\n[Use `summary(your.model.object)` function for more details.]\n")
}


#' @export
logLik.bizicount = function(object, ...) object$loglik

# AIC precomputed using k=2 in regression fn, so adjust it here according to AIC's k arg
#' @export
AIC.bizicount = function(object, ..., k=2) k/2*object$aic

#' @export
BIC.bizicount = function(object, ...) object$bic

#' @export
nobs.bizicount = function(object, ...) object$nobs

#' @export
vcov.bizicount = function(object, ...) object$covmat

#' @export
coef.bizicount = function(object, id=T, ...){
  if(id)
    return(object$coef)
  else
    return(object$coef.nid)
}

#' @export
fitted.bizicount = function(object, ...){

  with(object, {

  if(is.null(object$model))
    stop("bizicount() model must be fit with `keep=TRUE` to get fitted values.", call.=F)

  n.zi = sum(grepl("zi", margins))
  X = model[["X"]]
  Z = model[["Z"]]
  offset.ct = model[["offset.ct"]]
  offset.zi = model[["offset.zi"]]


  fit = list()

  for(i in 1:2){

    ct = coef[grepl(paste0("ct", i), names(coef))]

    if(grepl("zi", margins[i])){
      zi = coef[grepl(paste0("zi", i), names(coef))]
      pz = invlink.zi[[i]](Z[[i]] %*% zi + offset.zi[[i]])
      pc = (1-pz)*invlink.ct[[i]](X[[i]]%*%ct + offset.ct[[i]])
      fit[[i]] = pc

    } else {
      fit[[i]] = invlink.ct[[i]](X[[i]] %*% ct + offset.ct[[i]])

    }
  }

  out = do.call(cbind, fit)
  colnames(out) = paste0("fit_", object$outcomes)
  return(out)
  })
}

#' @export
summary.bizicount = function(object, ...){
  if(!any(class(object)=="bizicount"))
    stop("Object must be of class `bizicount`.")

  class(object) = c("summary.bizicount", "bizicount")
  return(object)
}

#' @export
print.summary.bizicount = function(x, stars=T, ...){
  old = getOption("show.signif.stars")
  options(show.signif.stars = stars)
  on.exit(options(show.signif.stars = old),
          add=T)

  with(x, {
  width = .618*getOption("width")
  zi = grep("zi", margins)
  nb = grep("nb", margins)
  n.zi = length(zi)
  n.nb = length(nb)

  pad = vector()

  # pad for printing of the coefmat for dispersion parameters
  if (n.nb > 0) {
    for (i in 1:2) {
      if (any(nb == i)){
        nds = rownames(coefmats[[paste0("disp_", outcomes[i])]])
        nct = rownames(coefmats[[i]])
        pad[i] = max(nchar(nct)) - max(nchar(nds)) + 1
        pad[pad < 0] = 0


        rownames(coefmats[[paste0("disp_", outcomes[i])]]) = paste0(nds, paste0(rep(" ", pad[i]), collapse=""))
      }
    }
  }

  #print call
  cat("Call:\n")
  cat(deparse(x$call), sep="\n")
  divider("=", width, prepend=T)

  for(i in 1:2){

    #count models
    modlabel(paste("Count Model:", outcomes[i]), append=T)
    printCoefmat(rbind(coefmats[[paste0("ct_", outcomes[i])]],
                       if(!is.null(coefmats[[paste0("disp_", outcomes[i])]]))
                         coefmats[[paste0("disp_", outcomes[i])]]
                       else
                         NULL
                       ),
                 signif.legend = if(i==2 && !any(zi==2) && stars) T else F)

    #zero inflated
    if ( any(zi == i) ){
      divider("+", .5*width, prepend=T)
      modlabel(paste("Zero Inflation:", outcomes[i]), append=T)
      printCoefmat(coefmats[[paste0("zi_", outcomes[i])]],
                   signif.legend = if(i == 2 && stars) T else F)
    }

    # dependence
    if ( i == 1 ){
      divider("-", width, prepend=T)
      printCoefmat(coefmats[["dependence"]],
                   signif.stars = if(stars) T else F,
                   signif.legend = F)
      divider("-", width, append=T)
    }
    else
      divider("=", width, prepend=T)
  }

  return(invisible(NULL))
  })

}


# Function for getting texreg output from bizicount objects.
#' @export
extract.bizicount = function(model, CI=NULL, id=T){
  if(!is.null(CI) && (CI > 1 || CI < 0) ) stop("`CI` must be between 0 and 1. ")
  if(!is.logical(id)) stop("`id` must be logical (T/F) value.")

  # Extract names of coefficients for each margin, create vector of new names that
  # that will match with output from zicreg so that new rows are not created in texreg output
  name = names(model$coef)

  pat = c("ct1_|zi1_",
          "ct2_|zi2_")
  nb = grepl("nb", model$margins)

  name.list = lapply(pat, function(x) grep(x, name, value=T))
  newname.list = lapply(name.list, function(x) gsub("^(.{2})[1-2](.*)", "\\1\\2", x))

  if(!id)
    newname.list = lapply(name.list, function(x) substr(x, 5, nchar(x)))

  for(i in 1:2){
    name.list[[i]] = c(name.list[[i]],
                       rep(paste0("disp.", i), nb[i]),
                       "dependence")
    newname.list[[i]] = c(newname.list[[i]],
                          rep("disp", nb[i]),
                          "dependence"
                          )
  }

  n = nobs(model)
  bic = BIC(model)
  aic = AIC(model)
  ll = logLik(model)

  gof = c(n, bic, aic, ll)

  gof.names = c("N", "BIC", "AIC", "LogLik")

  # create texreg output
  out = list()

  for(i in 1:2){
    out[[i]] = createTexreg(
      coef.names = newname.list[[i]],
      coef = model$coef[name.list[[i]]],
      se = model$se[name.list[[i]]],
      pvalues = model$p[name.list[[i]]],
      ci.low = if (!is.null(CI))
        model$coef[name.list[[i]]] - qnorm((1 + CI) / 2) * model$se[name.list[[i]]]
      else
        numeric(0),
      ci.up = if(!is.null(CI))
        model$coef[name.list[[i]]] + qnorm((1+CI)/2) * model$se[name.list[[i]]]
      else
        numeric(0),
      gof.names = gof.names,
      gof = gof,
      gof.decimal = c(F,T,T,T),
      model.name = paste0("Biv: ", model$outcomes[i])
    )
  }

  return(out)

}

#' @title The bizicount S4 Class
#' @description Note that `bizicount` objects are, in general, S3. However,
#' this S4 class is defined for compatability with \code{\link[texreg]{texreg}}.
#' Interaction with `bizicount` objects should generally use S3 syntax.
#' @slot coef Coefficients of the model
#' @slot coef.nid Coefficients without margin IDs
#' @slot coef.orig Coefficients prior to transformations, for Gaussian
#'   dependence and negative binomial dispersion.
#' @slot coef.orig.nid Coefficients prior to transforms, no margin IDs.
#' @slot se Asymptotic standard errors based on observed Fisher Information
#' @slot se.nid Standard errors without margin IDs
#' @slot z z-scores for parameter estimates
#' @slot z.nid z-scores without margin IDs
#' @slot p p-values for parameter estimates
#' @slot p.nid p-values without margin IDs
#' @slot coefmats A list containing coeficient matrices for each margin
#' @slot loglik Scalar log-likelihood at convergence
#' @slot grad Gradient vector at convergence
#' @slot n.iter Number of quasi-newton fitting iterations.
#' @slot covmat Covariance matrix of parameter estimates based on observed Fisher Information
#' @slot aic Model's Akaike information
#' @slot bic Model's Bayesian information criterion
#' @slot nobs Number of observations
#' @slot margins Marginal distributions used in fitting
#' @slot link.zi,link.ct Names of link functions used in fitting
#' @slot invlink.ct,invlink.zi Inverse link functions used in fitting (the
#'   actual function, not their names)
#' @slot outcomes Name of the response vector
#' @slot conv Integer telling convergence status.
#' @slot cop The copula used in fitting
#' @slot starts list of starting values used
#' @slot call The model's call
#' @slot model List containing model matrices, or `NULL` if `keep = F`.
#' @export
setClass("bizicount",
         representation(coef = "numeric",
                        coef.nid = "numeric",
                        coef.orig = "numeric",
                        coef.orig.nid = "numeric",
                        se = "numeric",
                        se.nid = "numeric",
                        z = "numeric",
                        z.nid = "numeric",
                        p = "numeric",
                        p.nid = "numeric",
                        coefmats = "numeric",
                        loglik = "numeric",
                        grad = "numeric",
                        n.iter = "integer",
                        covmat = "numeric",
                        aic = "integer",
                        bic = "numeric",
                        nobs = "integer",
                        margins = "character",
                        link.zi = "character",
                        link.ct = "character",
                        invlink.ct = "function",
                        invlink.zi = "function",
                        outcomes = "character",
                        conv = "integer",
                        cop = "character",
                        starts = "list",
                        call = "language",
                        model = "ANY"
         )
)

setMethod(texreg::extract,
          signature = "bizicount",
          definition= extract.bizicount)

#' @export
simulate.bizicount = function(object, nsim=250, seed=123, ...){
  if(is.null(object$model))
    stop("Must set `keep=T` in bizicount() to do diagnostics on model object.")

  if(!is.null(seed))
    set.seed(seed)

  with(object, {

  #observed = model$model[["y"]]
  n.zi = sum(grepl("zi", margins))
  nob = object$nobs
  X = model[["X"]]
  Z = model[["Z"]]
  offset.ct = model[["offset.ct"]]
  offset.zi = model[["offset.zi"]]
  disp1 = if(!is.null(coef["disp.1"])) coef["disp.1"] else NULL
  disp2 = if(!is.null(coef["disp.2"])) coef["disp.2"] else NULL

  psi = list()
  lam = list()
  disp = list()

  for(i in 1:2){

    ct = coef[grepl(paste0("^ct", i), names(coef))]
    lam[[i]] = as.vector(invlink.ct[[i]](X[[i]]%*%ct + offset.ct[[i]]))

    if(grepl("zi", margins[i])){
      zi = coef[grepl(paste0("^zi", i), names(coef))]
      psi[[i]] = as.vector(invlink.zi[[i]](Z[[i]] %*% zi + offset.zi[[i]]))
    } else
      NULL

    if(grepl("nb", margins[i]))
      disp[[i]] = coef[grepl(paste0("disp.", i), names(coef))]

  }

  margin.args = list()

  for (i in c(1, 2)) {
    margin.args[[i]] = switch(
      margins[i],
      "pois" = list(
        n      = nob,
        lambda = lam[[i]]
      ),
      "nbinom" = list(
        n     = nob,
        mu    = lam[[i]],
        size  = disp[[i]]
      ),
      "zip" = list(
        n      = nob,
        lambda = lam[[i]],
        psi    = psi[[i]]
      ),
      "zinb" = list(
        n     = nob,
        mu    = lam[[i]],
        size  = disp[[i]],
        psi   = psi[[i]]
      )
    )
  }

  sims = list()
  for(j in 1:2){
    sims[[j]] = replicate(nsim, do.call(paste0("r", margins[j]), margin.args[[j]]))

  }

  names(sims) = model$outcomes
  return(sims)

  }) # close out `with` command

}

#' @export
make_DHARMa = function(object, nsim=250, seed=123, method="PIT"){
  if( !any(class(object) == "bizicount") )
    stop("Function must be applied to bizicount class object.")

  sims = simulate(object, nsim=nsim, seed=seed)
  fit = fitted(object)

  dharmas = list()

  for(i in 1:2){
    dharmas[[i]] = createDHARMa(
      simulatedResponse = sims[[i]],
      observedResponse = object$model[["y"]][,i],
      fittedPredictedResponse = fit[,i],
      integerResponse = T,
      seed=FALSE,
      method=method
    )
  }

  return(dharmas)
}


## Generics for univariate regression model (zic.reg)
#' @export
summary.zicreg = function(object, ...) {
  if (!any(class(object) == "zicreg"))
    stop("Object must be of class `zicreg`.")

  class(object) = c("summary.zicreg", "zicreg")
  return(object)
}

# Print method for zicmod class object
#' @export
print.summary.zicreg = function(x, ...) {
  cat(paste("Call:", deparse1(x$call)), "\n")
  cat("\n", x$obj, paste0("(dist = ", x$dist, ")"), "\n")
  cat(
    "=========================================================================\n"
  )
  cat("Count Model", paste0("(link = ", x$link.ct, ")"),   "|\n")
  cat("--------------------------\n")
  printCoefmat(x$coefmat.ct, signif.legend = F)
  if (x$dist == "nbinom")
    cat("---\n\n")
  if (x$dist == "nbinom")
    printCoefmat(x$theta, signif.legend = F, dig.tst = 5)
  cat("-----------------------------------------------------------\n")
  cat("Zero Inflation", paste0("(link = ", x$link.zi, ")"), "|\n")
  cat("------------------------------\n")
  printCoefmat(x$coefmat.zi)
  cat(
    "=========================================================================\n"
  )
  cat("AIC: ", x$aic, "   BIC: ", x$bic, "\n")
}

#' @export
print.zicreg = function(x, ...) {
  cat(deparse1(x$call), "\n\n")
  cat("Count Model Coefficients:\n")
  print(x$coefmat.ct[, 1])
  cat("\nZero Inflation Coefficients:\n")
  print(x$coefmat.zi[, 1])

  if(!is.null(x$theta)){
    cat("\n")
    print(c(Theta=x$theta[1,1]))
  }

  cat("\n[Use `summary(your.model.object)` function for more details.]\n")

}


#' @export
coef.zicreg = function(object, ...) {
  n.ct = paste0("ct_", rownames(object$coefmat.ct))
  n.zi = paste0("zi_", rownames(object$coefmat.zi))
  cs   = c(object$coefmat.ct[, 1], object$coefmat.zi[, 1])
  names(cs) = c(n.ct, n.zi)
  cs
}

#' @export
logLik.zicreg = function(object, ...)
  object$loglik

#' @export
nobs.zicreg = function(object, ...)
  object$nobs

#' @export
BIC.zicreg = function(object, ...)
  object$bic

#' @export
AIC.zicreg = function(object, ..., k=2)
  2/k * object$aic

#' @export
vcov.zicreg = function(object, ...)
  object$covmat

#' @export
#' @name simulate.zicreg
#' @title This is a title
#' @description This is a description
simulate.zicreg = function(object, nsim = 250, seed=123, ...) {
  if (is.null(object$model))
    stop("Must set `keep=T` in zic.reg() to do diagnostics on model object.")

  if(!is.null(seed))
    set.seed(seed)

  with(object, {
    #observed = model$model[["y"]]
    nob = object$nobs
    X = model[["X"]]
    Z = model[["z"]]
    offset.ct = model[["offset.ct"]]
    offset.zi = model[["offset.zi"]]

    ct = as.vector(coefmat.ct[, 1])
    zi = as.vector(coefmat.zi[, 1])
    disp = if (!is.null(theta))
      theta[1]

    psi = make.link(link.zi)$linkinv(Z %*% zi + offset.zi)
    lam = make.link(link.ct)$linkinv(X %*% ct + offset.ct)

    dist.args = switch(
      dist,
      "pois" = list(
        n = nob,
        lambda = lam,
        psi = psi
      ),
      "nbinom" = list(
        n = nob,
        mu = lam,
        psi = psi,
        size = disp
      )
    )

    set.seed(seed)
    sims = replicate(nsim, do.call(if (grepl("p", dist))
      "rzip"
      else
        "rzinb", dist.args))

    return(sims)

  }) # close out `with` command

}

#' @export
fitted.zicreg = function(object, ...) {
  if (is.null(object$model))
    stop("`Keep` must be set to TRUE in original model call in order to obtain fitted mean.")

  with(object, {
    X = model[["X"]]
    Z = model[["z"]]
    offset.ct = model[["offset.ct"]]
    offset.zi = model[["offset.zi"]]

    ct = as.vector(coefmat.ct[, 1])
    zi = as.vector(coefmat.zi[, 1])

    psi = make.link(link.zi)$linkinv(Z %*% zi + offset.zi)
    lam = make.link(link.ct)$linkinv(X %*% ct + offset.ct)

    return(as.vector((1 - psi) * lam))
  })
}

#' @export
predict.zicreg = function(object,
                          newdata = NULL,
                          y.new = NULL,
                          type = "mean",
                          ...
) {

  type = match_arg(type, choices = c("mean", "prob", "psi", "lambda"))

  if (!is.null(newdata) && is.null(y.new) && type == "prob")
    stop("Must include arg `y.new` if `newdata` present and type='prob'.")

  if (nrow(newdata) != nrow(cbind(y.new)) && !is.null(y.new) && !is.null(newdata))
    stop("`newdata` and `y.new` not conformable.")

  if (is.null(newdata))
    warning(
      "`newdata` not specified; computations done on
            observed data (from original model)."
    )

  if (is.null(newdata) && type == "mean")
    return(fitted(object))

  with(object, {
    Z = model$z
    X = model$X
    ct = coefmat.ct[, 1]
    zi = coefmat.zi[, 1]
    offset.zi = if(!is.null(model$offset.zi)) model$offset.zi else 0
    offset.ct = if(!is.null(model$offset.ct)) model$offset.ct else 0

    disp = if (!is.null(theta))
      theta[1]

    if (is.null(newdata)) {
      psi.fit = make.link(link.zi)$linkinv(Z %*% zi + offset.zi)
      lam.fit = make.link(link.ct)$linkinv(X %*% ct + offset.ct)


      out = switch(type,
                   "prob" = switch(
                     dist,
                     "pois" = dzip(y, lam.fit, psi.fit),
                     "nbinom" = dzinb(
                       y,
                       size = disp,
                       mu = lam.fit,
                       psi = psi.fit
                     )
                   ),
                   "psi" = psi.fit,
                   "lambda" = lam.fit) #mean is returned at beginning of fn if newdata=NULL


    } else {

      coef.names = grep("Intercept",
                        names(coef(object)),
                        value = T,
                        invert = T)
      coef.names = sub("^.{3}", "", coef.names)
      assert(
        is.data.frame(newdata) &&   all(coef.names %in% names(newdata)),
        "`newdata` must be a dataframe containing all variables appearing in model (excluding intercept)."
      )

      X.new = newdata[, names(newdata) %in% names(ct)]
      Z.new = newdata[, names(newdata) %in% names(zi)]
      if (any(grepl("Intercept", names(ct))))
        X.new = cbind(1, X.new)
      if (any(grepl("Intercept", names(zi))))
        Z.new = cbind(1, Z.new)

      psi.new = make.link(link.zi)$linkinv(Z.new %*% zi + offset.zi)
      lam.new = make.link(link.ct)$linkinv(X.new %*% ct + offset.ct)

      out = switch(
        type,
        "mean" = (1 - psi.new) * lam.new,
        "prob" = switch(
          dist,
          "pois" = dzip(y.new, lam.new, psi.new),
          "nbinom" = dzinb(
            y.new,
            size = disp,
            mu = lam.new,
            psi = psi.new
          )
        ),
        "psi" = psi.new,
        "lambda" = lam.new
      )

    }

    return(out)

  })

}

# Function for texregging output
#' @name extract.zicreg
#' @title Texreg for zicreg objects
#' @description This is a method for the \code{\link[texreg]{extract}} generic
#'   to be used with \code{\link{zicreg-class}} objects that are output from the
#'   \code{\link{zic.reg}} function.
#' @method extract zicreg
#' @param model A zicreg model object (S3).
#' @param CI The two-tailed confidence level, if desired in the texreg object.
#' @param id Logical indicating whether to prepend equation identifiers to
#'   coefficient names (`ct_` for count parameters, `zi_` for zero-inflated parameters)
#' @return A \code{link[texreg]{texreg-class}} object, as produced by
#'   \code{link[texreg]{createTexreg}}, which can interface with all of that
#'   package's methods.
#' @author John Niehaus
#' @seealso \code{\link[texreg]{extract}}, \code{\link[texreg]{createTexreg}},
#'   \code{\link[bizicount]{zic.reg}}
#' @export
extract.zicreg = function(model, CI = NULL, id = T) {
  if (!is.null(CI) &&
      (CI > 1 || CI < 0))
    stop("`CI` must be between 0 and 1. ")
  if (!is.logical(id))
    stop("`id` must be a logical (T/F) value.")

  name = c(paste0("ct_", rownames(model$coefmat.ct)),
           paste0("zi_", rownames(model$coefmat.zi)),
           rownames(model$theta))

  if (!id)
    name = c(rownames(model$coefmat.ct),
             rownames(model$coefmat.zi),
             rownames(model$theta))

  colgrab = function(colnum) {
    c(model$coefmat.ct[, colnum],
      model$coefmat.zi[, colnum],
      model$theta[, colnum])
  }

  b = colgrab(1)
  se = colgrab(2)
  z  = colgrab(3)
  p = colgrab(4)
  n = nobs(model)
  bic = BIC(model)
  aic=AIC(model)
  ll = logLik(model)

  CI.low = if (!is.null(CI)) {
    b - qnorm((1 + CI) / 2) * se
  } else
    (numeric(0))

  CI.hi = if (!is.null(CI)) {
    b + qnorm((1 + CI) / 2) * se
  } else
    numeric(0)


  gof = c(n, bic, aic, ll)

  gof.names = c("N", "BIC", "AIC", "LogLik")


  tr = createTexreg(
    coef.names = name,
    coef = b,
    se = se,
    pvalues = p,
    gof.names = gof.names,
    gof = gof,
    ci.low = CI.low,
    ci.up = CI.hi,
    model.name = "ZIC Regression"
  )

  return(tr)
}

#' @title The zicreg S4 Class
#' @description Note that `zicreg` objects are, in general, S3. However,
#' this S4 class is defined for compatability with \code{\link[texreg]{texreg}}.
#' Interaction with `zicreg` objects should generally use S3 syntax, but the below
#' objects have the same name in both the S3 and S4 objects (but are in a list for S3).
#' @slot call The original function call
#' @slot obj The class of the object
#' @slot coef Vector of coefficients, with count, then zi, then dispersion.
#' @slot se Vector of asymptotic standard errors
#' @slot grad Gradient vector at convergence
#' @slot link.ct Name of link used for count portion
#' @slot link.zi Name of link used for zero-inflated portion
#' @slot dist Name of distribution used for count portion
#' @slot optimizer Name of optimization package used in fitting
#' @slot coefmat.ct Coefficient matrix for count portion
#' @slot coefmat.zi Coefficient matrix for zero-inflated portion
#' @slot coefmat.all Coefficient matrix for both parts of the model
#' @slot theta Coefficient matrix for dispersion, if applicable.
#' @slot covmat Asymptotic covariance matrix
#' @slot nobs Number of observations
#' @slot aic Akaike information
#' @slot bic Bayes information
#' @slot loglik Log-likelihood at convergence
#' @slot model List containing model matrices if `keep = TRUE`
#' @export
setClass(
  "zicreg",
  representation(
    call = "language",
    obj = "character",
    coef = "numeric",
    se = "numeric",
    grad = "ANY",
    link.ct = "character",
    link.zi = "character",
    dist = "character",
    optimizer = "character",
    coefmat.ct = "numeric",
    coefmat.zi = "numeric",
    coefmat.all = "numeric",
    theta = "ANY",
    covmat = "numeric",
    nobs = "integer",
    aic = "numeric",
    bic = "numeric",
    loglik = "numeric",
    convergence = "integer",
    model = "ANY"
  )
)
setMethod("extract", signature = "zicreg", definition = extract.zicreg)



