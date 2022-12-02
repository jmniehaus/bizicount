#' @name zi_test
#' @title He's (2019) test for zero-modification
#' @description This is an implementation of He et al. (2019)'s test for
#' zero-modification (discussed further in Tang & Tang (2019)). This is a test of
#' zero-*modification* instead of *inflation*, because the test is capable of detecting
#' both excessive or lack of zeros, but cannot determine the cause. For example, a mixed
#' data generating process could be generating structural zeros, implying a
#' zero-inflated distribution. However, overdispersion via a negative binomial
#' may also result in excessive zeros. Thus, the test merely determines whether
#' there are excessive (or lacking) zeros, but does not determine the process
#' generating this pattern. That in mind, typical tests in the literature are
#' inappropriate for zero-modified regression models, namely the Vuong, Wald,
#' score, and likelihood ratio tests. See the references below for more information
#' on this claim.
#'
#' @details The test compares the
#' observed proportion of zeros in the data to the expected proportion of zeros
#' under the null hypothesis of a Poisson distribution. This is done using
#' estimating equations to account for the fact that the expected proportion is
#' based on an estimated parameter vector, rather than the true parameter vector.
#' The test statistic is
#'
#' \eqn{\hat s = 1/n\sum_i (r_i - \hat p_i)}
#'
#' where \eqn{r_i = 1} if \eqn{y_i = 0}, otherwise \eqn{r_i = 0}, and \eqn{\hat p = dpois(0, exp(X\hat\beta)) = \hat E(r_i)}
#' is the estimated proportion of zeros under the assumption of a Poisson distribution
#' generated with covariates \eqn{X} and parameter vector \eqn{\hat\beta}.
#'
#' By the central limit theorem, \eqn{\hat s \sim AN(0, \sigma^2_s)}. However,
#' estimating \eqn{\hat \sigma_s} by a plug-in estimate using \eqn{\hat\beta} is inefficient
#' due to \eqn{\hat \beta} being an random variable with its own variance. Thus,
#' \eqn{\hat\sigma} is estimated via estimating equations in order to account for the
#' variance in \eqn{\hat \beta}.
#'
#' See the references below for more discussion and proofs.
#'
#' @param model A model object of class \code{\link{bizicount}} or \code{\link{glm}}.
#' If a \code{bizicount} model, then at least one margin must be specified as \code{"pois"}.
#' If a \code{glm} model, then the \code{family} must be \code{\link{poisson}}.
#' @param alternative The alternative hypothesis. One of \code{c("inflated", "deflated", "both")}.
#' These correspond to an upper tail, lower tail, or two-tailed test, respectively.
#' Default is \code{"inflated"}. Partial matching supported.
#' @references He, H., Zhang, H., Ye, P., & Tang, W. (2019). A test of inflated
#' zeros for Poisson regression models. Statistical methods in medical research,
#' 28(4), 1157-1169.
#'
#' Tang, Y., & Tang, W. (2019). Testing modified zeros for Poisson regression
#' models. Statistical Methods in Medical Research, 28(10-11), 3123-3141.
#' @author John Niehaus
#' @example inst/examples/zi_test_ex.R
#' @export
zi_test = function(model, alternative = 'inflated'){
     UseMethod("zi_test", model)
}

#' @export
zi_test.glm = function(model, alternative = 'inflated'){
     alternative = match_arg(alternative, choices = c("inflated", "deflated", "both"))
     assert(
          model$family$family == 'poisson',
          "glm family must be poisson()."
     )

     y = model.response(model.frame(model))
     y0 = (y == 0)
     V = vcov(model)
     lam = fitted(model, type = "response")
     d0 = dpois(0, lam)
     X = model.matrix(model)
     puX = (d0 * lam) %*% X

     tstat = sum(y0 - d0)/sqrt((sum(d0*(1-d0)) - puX %*% tcrossprod(V, puX)))

     pval = switch(
          alternative,
          "inflated" = pnorm(tstat, lower.tail = F),
          "deflated" = pnorm(tstat, lower.tail = T),
          "both"     = 2*pnorm(abs(tstat), lower.tail = F)
     )

     out = list(data.frame(
          H_a = alternative,
          Z_score = tstat,
          p_value = pval,
          n = length(y)
     ))

     rownames(out[[1]]) = colnames(model.frame(model))[1]
     class(out) = "zi_test"

     return(out)
}

#' @export
zi_test.bizicount = function(model, alternative = 'inflated'){
     alternative = match_arg(alternative, choices = c("inflated", "deflated", "both"))
     assert(
          any(model$margins == 'pois'),
          'At least one margin in bizicount model must be "pois".'
     )

     m_num = c(which(model$margins == "pois"))

     y = model$model[["y"]]

     y0 = (y == 0)
     lam = fitted(model)
     d0 = dpois(0, lam)
     X = model$model$X

     out = list()
     for(i in m_num){
          y0i = y0[,i]
          d0i = d0[,i]

          V = solve(crossprod(X[[i]]*lam[,i], X[[i]]))
          puX = (d0i * lam[,i]) %*% X[[i]]

          tstat = sum(y0i - d0i)/sqrt((sum(d0i*(1-d0i)) - puX %*% tcrossprod(V, puX)))
          pval = switch(
               alternative,
               "inflated" = pnorm(tstat, lower.tail = F),
               "deflated" = pnorm(tstat, lower.tail = T),
               "both"     = 2*pnorm(abs(tstat), lower.tail = F)
          )

          out[[i]] = data.frame(
               H_a = alternative,
               Z_score = tstat,
               p_value = pval,
               n = nrow(y)
          )
          rownames(out[[i]]) = colnames(y)[i]

     }

     class(out) = 'zi_test'
     return(out)
}

#' @export
print.zi_test = function(x, ...){
     cat("\n====================================================\n")
     cat("He et al. (2019)'s Test for Zero Modification\n")
     cat("--------------------------------------\n")
     cat("\nH_0:  Pr(y = 0 | x) = dpois(0 | x) \n\n\n")
     print(do.call(rbind.data.frame, x))
     cat("\n====================================================\n")
}
