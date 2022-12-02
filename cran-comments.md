# Bizicount Version 1.3.0 Comments

This is a minor release in the 1.*.0 series. Major changes include:

* Deprecated `scaling` parameter, as there is no reliable way to scale 
covariates properly when transformations and interactionsare introduced to model 
formulas.

* Deprecated `na.action` as there are no methods for `bizicount` that can
leverage alternative `na.actions` other than `na.omit`. 

See NEWS for other small details.

# R CMD check results

0 errors | 0 warnings | 1 note 

The release version flags a bad DOI: 10.1177/0962280217749991
I've manually checked and the DOI does work. Please advise if changes necessary.

# Test Environments 

* Ubuntu oldrel, rel, devel
* Windows oldrel, rel, devel
* R-winbuilder oldrel, rel, devel

# Revdep results 
`devtools::revdep` indicates no reverse dependencies.

