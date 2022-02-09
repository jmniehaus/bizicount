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
coef.bizicount = function(object, id = TRUE, ...){
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
print.summary.bizicount = function(x, stars = TRUE, ...){
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
#' @name extract.bizicount
#' @title Texreg for bizicount objects
#' @description This is a method for the \code{\link[texreg]{extract}} generic
#'   to be used with objects that are output from the \code{\link{bizicount}}
#'   function. The results can be used with any of the
#'   \code{\link[texreg]{texreg-package}} generics.
#' @method extract bizicount
#' @param model A \code{\link{bizicount-class}} model object (S3).
#' @param CI The two-tailed confidence level, if confidence intervals are
#'   desired in the texreg object, otherwise `NULL`.
#' @param id Logical indicating whether to prepend equation identifiers to
#'   coefficient names (`ct_` for count parameters, `zi_` for zero-inflated parameters)
#' @return A \code{\link[texreg]{texreg-class}} object, as produced by
#'   \code{\link[texreg]{createTexreg}}, which can interface with all of that
#'   package's generics.
#' @note Users can typically just call \code{\link[texreg]{texreg}} directly on
#'   a \code{\link{bizicount-class}} object, instead of first extracting and
#'   then calling texreg.
#' @example inst/examples/extract_bizicount_ex.R
#'
#' @references Leifeld, Philip (2013). texreg: Conversion of Statistical Model
#'   Output in R to LaTeX and HTML Tables. Journal of Statistical Software,
#'   55(8), 1-24. URL http://dx.doi.org/10.18637/jss.v055.i08.
#'
#' @author John Niehaus
#' @seealso \code{\link[texreg]{extract}}, \code{\link[texreg]{createTexreg}},
#'   \code{\link[bizicount]{bizicount}}
#' @export
extract.bizicount = function(model, CI = NULL, id = TRUE){
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
#' @description Note that `bizicount` objects are generally S3, and should use
#'   S3 syntax. This S4 class is defined only for compatability with
#'   \code{\link[texreg]{texreg}}. However, the contents of `bizicount` objects
#'   is the same in both S3 and S4, so the descriptions below apply in both cases.
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


#' @name simulate.bizicount
#' @title Simulating response values using parameters from fitted bizicount models
#' @description Simulates random response values using the fitted conditional
#'   mean function for each margin of a \code{\link{bizicount-class}} object.
#'   Primarily for use with the \code{\link[DHARMa]{DHARMa}} package.
#' @param object A fitted \code{\link{bizicount-class}} object, as returned by
#'   \code{\link{bizicount}}.
#' @param nsim Number of simulated response values from the fitted model. E.g.,
#'   `nsim = 250` will simulate each observation 250 times, for \eqn{n \times
#'   250} total observations.
#' @param seed Seed used for simulating from fitted model. If `NULL`, no seed is
#'   set.
#' @param ... Ignored.
#' @return A length 2 list, with each entry containing a numeric \eqn{n X
#'   nsim} matrix for each margin of the bizicount model. Rows index
#'   the observation, and columns index the simulated dataset number.
#' @example inst/examples/simulate_bizicount_ex.R
#'
#' @references Florian Hartig (2022). DHARMa: Residual Diagnostics for
#'   Hierarchical (Multi-Level / Mixed) Regression Models. R package version
#'   0.4.5. https://CRAN.R-project.org/package=DHARMa
#' @author John Niehaus
#' @seealso \code{\link[DHARMa]{createDHARMa}}, \code{\link[DHARMa]{simulateResiduals}}
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

#' @name make_DHARMa
#' @title DHARMa-class objects from bizicount models
#' @description A wrapper around the \code{\link[DHARMa]{DHARMa}} package's
#'   \code{\link[DHARMa]{createDHARMa}} function. Creates a list of DHARMa
#'   objects, one for each margin of a \code{\link{bizicount-class}} object, using
#'   simulated responses from \code{\link{simulate.bizicount}}.
#' @seealso \code{\link[DHARMa]{DHARMa}}, \code{\link[DHARMa]{createDHARMa}},
#'   \code{\link{simulate.bizicount}}
#' @author John Niehaus
#' @note This is merely a wrapper around the \code{\link[DHARMa]{createDHARMa}}
#'   function to conveniently get DHARMa objects for each margin of a bizicount
#'   model.
#' @return A list of \code{\link[DHARMa]{DHARMa}} objects.
#' @param object A \code{\link{bizicount-class}} object, as returned by \link{bizicount}.
#' @param nsim Number of simulated responses from the fitted model to use for diagnostics.
#' @param seed Random seed for simulating from fitted model.
#' @param method See \code{\link[DHARMa]{createDHARMa}}.
#' @example inst/examples/make_dharma_ex.R
#' @references Florian Hartig (2022). DHARMa: Residual Diagnostics for
#'   Hierarchical (Multi-Level / Mixed) Regression Models. R package version
#'   0.4.5. https://CRAN.R-project.org/package=DHARMa
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
