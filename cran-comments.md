# Resubmission
This is a resubmission. In this version, I have:

* Added supporting references with DOIs to DESCRIPTION. 

* Amended function definitions to use TRUE/FALSE rather than T/F.

* Added a return value to `predict.zicreg.Rd`.

* Added more detailed return values to `bizicount.Rd` and `zic.reg.Rd`.

* Removed `options()` call in `bizicount.R`, wrapped relevant calls in `suppressWarnings()` 
instead. 

Two additional requests were made, for which no changes were implemented.
Explanations include:

* It was requested that `\dontrun` be changed to `\donttest` in two examples 
(`dzip.Rd` and `dzinb.Rd`). The two examples are there to demonstrate that an 
error will occur when arguments lack a certain property. Thus, the test intentionally
throws an error. When running R CMD check, these examples flag an error when 
wrapped in `\donttest`, but do not throw the error when wrapped in `\dontrun`. 
Therefore, `\dontrun` has been kept in order to have complete examples for users
while still passing CMD check.

* The previous submission's reply indicated that I need to include 
all relevant copyright holders. For example, Diane Lambert is cited in the 
references section of `man/zic.reg.Rd`. I am the sole author and copyright holder
for all code found in this package. I merely cited to Diane Lambert's paper on
zero-inflated count models as a reference for the statistical theory behind
the method; no code, statements, or otherwise were used in any way from her work, 
or any other work. Therefore, I do not know of any changes that need to be made
to authors, contributors, or copyright holders, as all work here is mine alone.

# Test Environments 

* Local Ubuntu 20.04, R-release
* Local Windows 10 Pro, R-release
* Remote MacOS 11 runner, R-release
* Remote MacOS 11 runner, R-devel
* Remote Ubuntu 20.04, R-release
* Remote Ubuntu 20.04 R-devel
* Win-builder, R-release
* Win-builder, R-devel

# R CMD check results
0 errors | 0 warnings | 1 note 

This note occurs because this is a new package. 

# Revdep results 
This is the first release; there are no reverse dependencies.

