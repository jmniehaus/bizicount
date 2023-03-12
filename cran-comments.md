# Bizicount Version 1.3.1 Comments

This is a bugfix release to fix R-CMD-CHECK issues found here:
https://cran.r-project.org/web/checks/check_results_bizicount.html.

Two issues were raised:

     1. A warning about generic methods 
     
     2. An error stating that `rlang` could not be loaded. 
     
The issue about generic methods has been resolved in this release. 

The issue about `rlang` has not been reproduced on the indicated platform. I 
ran checks on R-winbuilder-devel and received no errors. Thus, I am assuming that
the namespace issue was server-side or due to an issue with a dependency that 
has since been resolved. 

# R CMD check results

0 errors | 0 warnings | 0 notes

# Test Environments 

* Ubuntu oldrel, rel, devel
* Windows oldrel, rel, devel
* R-winbuilder oldrel, rel, devel

# Revdep results 
`devtools::revdep` indicates no reverse dependencies.

