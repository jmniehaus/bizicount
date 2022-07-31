test_that("univariate dists work with correct params",{

     expect_equal(dzip(5,3,.5), .5*dpois(5,3))

     expect_equal(dzip(0, 5, 0.5), 0.5 + 0.5*dpois(0,5))

     expect_equal(qzip(.25, 5, .5), 0)

     expect_equal(qzip(.25, 5, .1), 3)

     expect_equal(pzip(5, 3, .5), .5 + .5*ppois(5, 3))

     expect_equal(pzip(Inf, 3, .5), 1)
     }
)

test_that("univariate dists work with incorrect params", {

     expect_warning(val <- dzip(-5, 5, .2), "Negative")
     expect_true(val == 0)

     expect_warning(dzip(1.2, 5, .2), "non-integer")

     expect_warning(val <- dzip(c(1,2,3), c(-5, 1, -2), c(.2, .5, .3)), "NaNs")
     expect_true(sum(is.na(val)) == 2)

     expect_warning(val <- dzip(c(1,2,3), c(5, 1, 2), c(1.2, .5, -.3)), "NaNs")
     expect_true(sum(is.na(val)) == 2)

     expect_warning(val <- dzip(c(1,2,3), c(-5, -1, 2), c(1.2, .5, -.3)), "NaNs")
     expect_true(all(is.na(val)))

     expect_warning(val <- dzip(c(1,-2,3), c(5, -1, 2), c(.2, .5, -.3)), "Negative|NaNs")
     expect_true(sum(is.na(val)) == 2)


     expect_warning(val <- pzip(-5, 5, .2), "Negative")
     expect_true(val == 0)

     expect_warning(val <- pzip(c(1,2,3), c(-5, 1, -2), c(.2, .5, .3)), "NaNs")
     expect_true(sum(is.na(val)) == 2)

     expect_warning(val <- pzip(c(1,2,3), c(5, 1, 2), c(1.2, .5, -.3)), "NaNs")
     expect_true(sum(is.na(val)) == 2)

     expect_warning(val <- pzip(c(1,2,3), c(-5, -1, 2), c(1.2, .5, -.3)), "NaNs")
     expect_true(all(is.na(val)))

     expect_warning(val <- pzip(c(1,-2,3), c(5, -1, 2), c(.2, .5, -.3)), "Negative|NaNs")
     expect_true(sum(is.na(val)) == 2)


     expect_warning(val <- qzip(c(.1, .2, .3), c(-5, 1, -2), c(.2, .5, .3)), "NaNs")
     expect_true(sum(is.na(val)) == 2)

     expect_warning(val <- qzip(c(.1, .2, .3), c(5, 1, 2), c(1.2, .5, -.3)), "NaNs")
     expect_true(sum(is.na(val)) == 2)

     expect_warning(val <- qzip(c(1.5, .5, -.2), 5, .2), "NaNs")
     expect_true(sum(is.na(val)) == 2)


     expect_warning( rzip(1, c(-5, 5), .2), "NaNs")

     expect_warning( rzip(3, 5, c(1.2, .2, -.3)), "NaNs")


     expect_warning(val <- rzinb(3, size = c(-1, .5, 1.2), mu = c(1,2,3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 1 )

     expect_warning(val <- rzinb(3, size = c(1, .5, 0), mu = c(-1,2,-3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- rzinb(3, size = c(1, -1, 0), mu = c(-1,2,-3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 3 )

     expect_warning(val <- rzinb(3, size = c(-1, 1, -1), mu = c(1,-2,3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 3 )

     expect_warning(val <- rzinb(3, size = c(1, 1, 1), mu = c(1,2,3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- rzinb(3, size = c(-1, 1, -1), mu = c(1,2,3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- rzinb(3, size = c(1, 1, 1), mu = c(-1,2,-3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- rzinb(3, size = c(-1, 1, -1), mu = c(-1,2,-3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)



     expect_warning(val <- dzinb(c(-1, 2, 3), size = c(1, .5, .1), mu = c(1,2,3), psi = c(.1,.2,.3)), "Negative")
     expect_warning(val <- dzinb(c(1.5, 2, 3), size = c(1, .5, .1), mu = c(1,2,3), "non-integer"))

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(-1, .5, 1.1), mu = c(1,2,3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 1 )

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(1, .5, 1), mu = c(-1,2,-3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(1, -1, 0), mu = c(-1,2,-3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 3 )

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(-1, 1, -1), mu = c(1,-2,3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 3 )

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(1, 1, 1), mu = c(1,2,3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(-1, 1, -1), mu = c(1,2,3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(1, 1, 1), mu = c(-1,2,-3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- dzinb(c(-1, 2, -3), size = c(-1, 1, -1), mu = c(-1,2,-3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)


     expect_warning(val <- pzinb(c(-1, 2, 3), size = c(1, .5, .1), mu = c(1,2,3), psi = c(.1,.2,.3)), "Negative")

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(-1, .5, 1.1), mu = c(1,2,3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 1 )

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(1, .5, 1), mu = c(-1,2,-3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(1, -1, 0), mu = c(-1,2,-3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 3 )

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(-1, 1, -1), mu = c(1,-2,3),  psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 3 )

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(1, 1, 1), mu = c(1,2,3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(-1, 1, -1), mu = c(1,2,3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(1, 1, 1), mu = c(-1,2,-3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)

     expect_warning(val <- pzinb(c(-1, 2, -3), size = c(-1, 1, -1), mu = c(-1,2,-3), psi = c(-.1,.2,-.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2)



     expect_warning(val <- qzinb(c(-.5, .2, 1.1), size = c(1, .5, .1), mu = c(1,2,3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- qzinb(c(.5, .2, .1), size = c(-1, .5, 1.1), mu = c(1,2,3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 1 )

     expect_warning(val <- qzinb(c(.5, .2, .1), size = c(1, .5, .1), mu = c(-1,2,-3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- qzinb(c(.5, .2, .1), size = c(1, .5, .1), mu = c(1,2,3), psi = c(-.1,.2,1.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- qzinb(c(-.5, .2, 1.1), size = c(-1, .5, 1.1), mu = c(1,2,3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- qzinb(c(-.5, .2, 1.1), size = c(1, .5, 1.1), mu = c(-1,2,-3), psi = c(.1,.2,.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- qzinb(c(-.5, .2, 1.1), size = c(1, .5, 1.1), mu = c(1,2,3), psi = c(-.1,.2,1.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- qzinb(c(-.5, .2, 1.1), size = c(-1, .5, 1.1), mu = c(1,2,3), psi = c(-.1,.2,1.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )

     expect_warning(val <- qzinb(c(-.5, .2, 1.1), size = c(-1, .5, 1.1), mu = c(-1,2,-3), psi = c(-.1,.2,1.3)), "NaNs")
     expect_true( sum(is.na(val)) == 2 )





     })


test_that("all bizicount margins work", {
     set.seed(123)
     y1 = rzinb(500, .3, .05, 8)
     y2 = rzinb(500, .15, .1, 5)

     grid = expand.grid(c("pois", "nbinom", "zip", "zinb"), c("pois", "nbinom", "zip", "zinb"), stringsAsFactors = F)
     grid$f1 = ifelse(grepl("zi", grid[,1]), "y1 ~ 1 | 1", "y1 ~ 1")
     grid$f2 = ifelse(grepl("zi", grid[,2]), "y2 ~ 1 | 1", "y2 ~ 1")


     out = list()
     for(i in seq_len(nrow(grid))){
          out[[i]] = suppressWarnings(tryCatch(
               bizicount(as.formula(grid$f1[i]),
                         as.formula(grid$f2[i]),
                         margins = unlist(grid[i, 1:2])),
               error = function(e)
                    NULL
          ))
     }

     expect_false(any(sapply(out, is.null)))

})
