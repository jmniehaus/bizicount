# Bizicount Version 1.3.4 Comments

* Fixed crossrefs per latest notes on CRAN build.

# R CMD check results

0 errors | 0 warnings | 0/1 notes (depending on test env)

## Note Discussion

On some builds, the JSTOR URL from the DESCRIPTION file returns a 403, but I've checked that this link is valid and it is. 
My guess is that the check process triggers some robots.txt blocking the request, but I can't be certain. Please advise if this should be pursued further.

# Test Environments 

* Ubuntu oldrel, rel, devel
* R-winbuilder oldrel, rel, devel

# Revdep results 
`devtools::revdep` indicates no reverse dependencies.

