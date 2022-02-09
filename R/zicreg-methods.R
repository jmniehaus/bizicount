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
#' @title Simulating response values from fitted univariate zero-inflated count
#'   regression model
#' @description Simulates responses using the fitted parameters from a
#'   \code{\link{zicreg-class}} object, as returned by \code{\link{zic.reg}}.
#'   Primarily useful for methods found in \code{\link[DHARMa]{DHARMa}} package. See 'Examples.'
#' @param object A \code{\link{zicreg-class}} omodel bject, as returned by \code{\link{zic.reg}}.
#' @param nsim Number of simulated datasets to create.
#' @param seed Random seed for random number generation in simulations. If
#'   `NULL`, no seed is set.
#' @param ... Ignored.
#' @returns A numeric \eqn{n x nsim} matrix, with rows indexing
#'   observations, and columns indexing the simulation number.
#' @example /inst/examples/simulate_zicreg_ex.R
#' @references Florian Hartig (2022). DHARMa: Residual Diagnostics for
#'   Hierarchical (Multi-Level / Mixed) Regression Models. R package version
#'   0.4.5. https://CRAN.R-project.org/package=DHARMa
#' @author John Niehaus
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


#' @name predict.zicreg
#' @title Predictions for univariate zero-inflated count regression models
#' @description Predicts the mean, probability, count mean, or zero-inflation
#'   probability for new data using parameters from a fitted zero-inflated count
#'   regression model.
#'
#' @param object A fitted \code{\link{zic.reg}} object.
#'
#' @param newdata A \code{\link[base]{data.frame}} containing new values of the same covariates appearing in fitted model.
#'
#' @param y.new An optional vector of new response values, used only for `type = "prob"`.
#'
#' @param type String, one of `c("mean", "prob", "psi", "lambda")`. `"mean"`
#'   will predict the conditional mean of the mixture distribution, `"prob"`
#'   will predict the probability of a new response value, `"psi"` will predict
#'   the probability of zero-inflation, and `"lambda"` will predict the mean of
#'   the count portion of the mixture distribution. NOTE: Setting `type =
#'   "mean"` and leaving `newdata = NULL` is the same as calling
#'   `fitted(object)`.
#'
#' @param ... Ignored.
#'
#' @return A numeric vector containing the predictions using the model parameters.
#'
#'
#'
#' @example /inst/examples/predict_zicreg_ex.R
#'
#' @author John Niehaus
#'
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
#'   to be used with objects that are output from the \code{\link{zic.reg}}
#'   function. The results can then interface with the
#'   \code{\link[texreg]{texreg-package}}, as shown in examples below.
#' @method extract zicreg
#' @param model A zicreg model object, returned by \code{\link{zic.reg}}.
#' @param CI The two-tailed confidence level, if desired in the resulting
#'   \code{\link[texreg]{texreg}} object.
#' @param id Logical indicating whether to prepend equation identifiers to
#'   coefficient names (`ct_` for count parameters, `zi_` for zero-inflated parameters)
#' @return A \code{\link[texreg]{texreg-class}} object, as produced by
#'   \code{\link[texreg]{createTexreg}}, which can interface with all of that
#'   package's generics. See 'Examples.'
#' @author John Niehaus
#' @seealso \code{\link[texreg]{extract}}, \code{\link[texreg]{createTexreg}},
#'   \code{\link[bizicount]{zic.reg}}
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software,
#'   55(8), 1-24. URL http://dx.doi.org/10.18637/jss.v055.i08.
#'
#'
#' @example /inst/examples/extract_zicreg_ex.R
#' @export
extract.zicreg = function(model, CI = NULL, id = TRUE) {
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
#' @slot convergence Convergence code from optimization routine.
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
setMethod(texreg::extract, signature = "zicreg", definition = extract.zicreg)



