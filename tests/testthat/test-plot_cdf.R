library("SimMultiCorrData")
context("Plots of theoretical cdf")

skip_on_cran()
W <- calc_theory("Weibull", c(3, 5))
Six <- seq(0.1, 1, 0.01)
CF <- find_constants("Fleishman", W[3], W[4])$constants
CP <- find_constants("Polynomial", W[3], W[4], W[5], W[6], Six = Six)$constants

test_that("works for Fleishman method", {
  expect_is(plot_cdf(CF, "Fleishman", mu = W[1], sigma = W[2],
    calc_cprob = FALSE), "ggplot")
  expect_is(plot_cdf(CF, "Fleishman", mu = W[1], sigma = W[2],
    calc_cprob = TRUE, delta = 5, lower = 0, upper = 1), "ggplot")
})

test_that("works for Polynomial method", {
  expect_is(plot_cdf(CP, "Polynomial", mu = W[1], sigma = W[2],
    calc_cprob = FALSE), "ggplot")
  expect_is(plot_cdf(CP, "Polynomial", mu = W[1], sigma = W[2],
    calc_cprob = TRUE, delta = 5, lower = 0, upper = 1), "ggplot")
})
