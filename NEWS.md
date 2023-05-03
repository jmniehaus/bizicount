# News

## bizicount 1.3.2

### Major 

* None

### Minor 

* Added terrorism data from associated journal article. Can be accessed by the 
object called `terror` which is lazily loaded upon loading the package.

## bizicount 1.3.1

### Major 

* None 

### Minor 

* There was a warning in R-CMD-CHECK on generic method consistency, which was 
fixed. 

## bizicount 1.3.0

* Deprecated `scaling` parameter, as there is no reliable way to scale 
covariates properly when transformations and interactionsare introduced to model 
formulas.

* Deprecated `na.action` as there are no methods for `bizicount` that can
leverage alternative `na.actions` other than `na.omit`. 

## bizicount 1.2.0 

### Major 

* Add `zi_test` function to implement tests for zero-modification as found in 
He et al. (2019).

* Fix bug in `scaling` parameter where interaction terms would be scaled, not
just their constituent covariates. 

* Fix bug causing standard error on dependence parameter in Gaussian copula 
to be incorrect. 


### Minor 

* Change default value of `keep` parameter in `bizicount()` function to `TRUE`
so that model matrices, etc., are stored in the output object by default. 

* Update documentation with new methods, typo fixes.

* Refactor, add more tests, and minor bug fixes for univariate distribution functions.




***

## bizicount 1.1.0
### Minor
* Fixed improper univariate distribution aliasing in documentation (#1)
* Added more badges to readme.md (#4)
* Pointed BugReports link to correct GitHub location (#5)
* Added new test to prevent #6 from happening again. 

### Major
* Fixed fatal error that occurred when the first marginal distribution was zero-inflated,
 while second margin was not (#6). 
* Fixed error when trying to use NumDeriv to get hessian if NLM hessian failed (#8).
* Added a parameter to `bizicount()` that scales continuous covariates automatically
for users (#7). 

***

## bizicount 1.0.0

The first release of the package, support for Gaussian and Frank copulas,
zero-inflated and non-inflated Poisson and negbin distributions. Extended
DHARMa and texreg generics for bizicount and zicreg objects.

