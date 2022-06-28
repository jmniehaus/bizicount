# News

# bizicount 1.0.0.9000 (dev)

## Minor
* Fixed improper univariate distribution aliasing in documentation (#1)
* Added more badges to readme.md (#4)
* Pointed BugReports link to correct GitHub location (#5)
* Added new test to prevent #6 from happening again. 

## Major
* Fixed fatal error that occurred when the first marginal distribution was zero-inflated,
 while second margin was not (#6). 
* Fixed error when trying to use NumDeriv to get hessian if NLM hessian failed (#8).
* Added a parameter to `bizicount()` that scales continuous covariates automatically
for users (#7). 

# bizicount 1.0.0

The first release of the package, support for Gaussian and Frank copulas,
zero-inflated and non-inflated Poisson and negbin distributions. Extended
DHARMa and texreg generics for bizicount and zicreg objects.

