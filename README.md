# bizicount

The `bizicount` R-package is primarily for estimating copula-based bivariate 
count regression models, with and without zero-inflation or overdispersion. 
However, it also includes functions to do the following:
     * Estimate univariate zero-inflated count models (`zic.reg`)
     * Evaluate univariate zero-inflated count distribution CDFs, PDFs, and 
     and quantile functions, and generate random zero-inflated counts 
     (`pzip` or `pzinb`, `dzip` or `dzinb`, `qzip` or `qzinb`, and `rzip` or `rzinb`)
     * Carry out post estimation diagnostics using the [`DHARMa`](https://github.com/florianhartig/DHARMa) package (`simulate.bizicount`, `simulate.zicreg`, `make_DHARMa`)
     * Produce professional tables in latex, word, or plain-text using the [`texreg`](https://github.com/leifeld/texreg) package (`extract.bizicount`, `extract.zicreg`). 
     
## Installation

To install from CRAN: 
`install.packages("bizicount")`

To install from GitHub:
`devtools::install_github("jmniehaus/bizicount", dependencies = TRUE)`

# Acknowledgements 
