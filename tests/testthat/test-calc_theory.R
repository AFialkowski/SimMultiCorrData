library("SimMultiCorrData")
context("Theoretical standardized cumulants")

skip_on_cran()
tol <- 1e-5

test_that("works given a distribution name", {
  expect_equal(all.equal(calc_theory(Dist = "Chisq", params = 2)[6], 120,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(calc_theory(Dist = "Laplace", params = c(0, 1))[6],
    30, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(calc_theory(Dist = "Dagum", params = c(1, 1, 2))[6],
    5.487257, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(calc_theory(Dist = "Beta-Normal",
    params = c(0.1, 4, 2, 1))[6], -1.2950460, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works given a pdf fx", {
  expect_equal(all.equal(calc_theory(fx = function(x) 0.5 * dbeta(x, 6, 3) +
    0.2 * dbeta(x, 4, 1.5) + 0.3 * dbeta(x, 10, 20), lower = 0, upper = 1)[6],
    6.5350069, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(calc_theory(fx = function(x) 0.4 * dnorm(x, -2, 1) +
    0.6 * dnorm(x, 2, 1), lower = -Inf, upper = Inf, sub = 500)[6],
    6.1732675, tolerance = tol, check.attributes = FALSE), TRUE)
})
