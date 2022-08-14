# Bizicount Version 1.2.0 Comments

This is a minor release in the 1.*.0 series. Major changes include:

* Adding a new hypothesis test for zero-modification for glm objects and 
bizicount objects. 

* Bug fix for incorrect standard error on the dependence parameter to the 
gaussian copula

* The internal scaling function was incorrectly scaling interaction terms, 
not just their constituent variables. 

See NEWS for other small details.

# R CMD check results

0 errors | 0 warnings | 0 note 

# Test Environments 

* Ubuntu oldrel, rel, devel
* Windows oldrel, rel, devel
* R-winbuilder devel

# Revdep results 
`devtools::revdep` indicates no reverse dependencies.

