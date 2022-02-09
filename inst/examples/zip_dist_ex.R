# Unequal lengths, but one of them is length 1, others are same length (3).
# No error.

x = c(1,2,3)
lambda = c(3,4,5)
psi = .1

dzip(x, lambda, psi)


# unequal lengths, at least one of them is not length 1,
# error

\dontrun{

x = c(1,2,3)
lambda = c(3,4)
psi = .1

dzip(x, lambda, psi)

}

# unequal lengths, at least one of them is not length 1.
# but set recycle = T to permit arbitrary recycling.

x = c(1,2,3)
lambda = c(3,4)
psi = .1

dzip(x, lambda, psi, recycle=TRUE)
