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

     expect_warning(val <- dzip(5, -5, .2), "NaNs")
     expect_true(is.nan(val))

     expect_warning(val <- dzip(5, 5, -.2), "NaNs")
     expect_true(is.nan(val))

     expect_warning(val <- pzip(-5, 5, .2), "Negative")
     expect_true(val == 0)

     expect_warning(val <- pzip(5,-5,.2), "NaNs")
     expect_true(is.nan(val))

     expect_warning(val <- pzip(5, 5, 5), "NaNs")
     expect_true(is.nan(val))

     expect_warning(val <- qzip(.5, -5, .2), "NaNs")
     expect_true(is.nan(val))

     expect_warning(val <- qzip(.5, 5, -.2), "NaNs")
     expect_true(is.nan(val))

     expect_warning(val <- qzip(1.5, 5, .2), "NaNs")
     expect_true(is.nan(val))

     expect_warning(val <- rzip(1, -5, .2), "NaNs")

     expect_warning(val <- rzip(1, 5, 1.2), "NaNs")


     }
)
