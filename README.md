# bizicount

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/jmniehaus/bizicount/workflows/R-CMD-check/badge.svg)](https://github.com/jmniehaus/bizicount/actions)
  [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bizicount?color=lightgrey)](https://cran.r-project.org/package=bizicount)
  [![minimal R version](https://img.shields.io/badge/R%3E%3D-4.1.0-6666ff.svg)](https://cran.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Downloads:](https://cranlogs.r-pkg.org/badges/grand-total/bizicount?color=red)](https://cran.r-project.org/package=bizicount)
[![Codecov test coverage](https://codecov.io/gh/jmniehaus/bizicount/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jmniehaus/bizicount?branch=main)
  <!-- badges: end -->


The `bizicount` R-package is primarily for estimating copula-based bivariate 
count regression models, with and without zero-inflation or overdispersion. However,
its full suite of functions can:

* Estimate copula-based bivariate zero-inflated count regression models, as well as 
non-inflated models. 

* Estimate univariate zero-inflated count models (`zic.reg`)
     
* Carry out post estimation diagnostics using the [`DHARMa`](https://github.com/florianhartig/DHARMa) package (`simulate.bizicount`, `simulate.zicreg`, `make_DHARMa`)
     
* Produce professional tables in latex, word, or plain-text using the [`texreg`](https://github.com/leifeld/texreg) package (`extract.bizicount`, `extract.zicreg`). 

* Test for zero modification using `zi_test()` [(He et al. 2019).](https://pmc.ncbi.nlm.nih.gov/articles/PMC6345607/)

* Evaluate univariate zero-inflated count distribution CDFs, PDFs, and quantile 
functions, and generate random zero-inflated counts (`pzip` or `pzinb`, `dzip` or `dzinb`, `qzip` or `qzinb`, and `rzip` or `rzinb`)
     
## Installation

To install from CRAN: 
```install.packages("bizicount")```

To install from GitHub:
```devtools::install_github("jmniehaus/bizicount", dependencies = TRUE)```

