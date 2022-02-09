# zero-inflated negative binomial examples

# two unique lengthed arguments, one is length 1 though. No error.

dzinb(4, size=.25, mu= c(1,2,3), psi=c(.2, .1, .15))


# two unique lengthed arguments, one of them is not length 1
# error
\dontrun{

     dzinb(5, size=c(.25, .3), mu= c(1,2,3), psi=c(.2, .1, .15))

}


# two unique lengthed arguments, one of them is not length 1, set
# recycle = T, no error but can give innacurate results.

dzinb(5, size=c(.25, .3), mu= c(1,2,3), psi=c(.2, .1, .15), recycle=TRUE)
