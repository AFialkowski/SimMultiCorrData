library("SimMultiCorrData")
context("Plots of simulated cdf")

skip_on_cran()
W <- calc_theory("Weibull", c(3, 5))
Six <- seq(0.1, 1, 0.01)
WF <- nonnormvar1(method = "Fleishman", means = W[1], vars = W[2]^2,
                  skews = W[3], skurts = W[4])$continuous_variable
WP <- nonnormvar1(method = "Polynomial", means = W[1], vars = W[2]^2,
                  skews = W[3], skurts = W[4], fifths = W[5],
                  sixths = W[6], Six = Six)$continuous_variable

marginal <- list(0.3)
lam <- 0.5
size <- 2
prob <- 0.75
Rey <- matrix(0.4, 3, 3)
diag(Rey) <- 1
Sim1 <- rcorrvar(k_cat = 1, k_pois = 1, k_nb = 1, marginal = marginal,
  support = list(c(0, 1)), lam = lam, size = size, prob = prob, rho = Rey)

test_that("works for Fleishman method", {
  expect_is(plot_sim_cdf(WF[, 1], calc_cprob = FALSE), "ggplot")
  expect_is(plot_sim_cdf(WF[, 1], calc_cprob = TRUE, delta = 5), "ggplot")
})

test_that("works for Polynomial method", {
  expect_is(plot_sim_cdf(WP[, 1], calc_cprob = FALSE), "ggplot")
  expect_is(plot_sim_cdf(WP[, 1], calc_cprob = TRUE, delta = 5), "ggplot")
})

test_that("works for ordinal variable", {
  expect_is(plot_sim_cdf(Sim1$ordinal_variables[, 1], calc_cprob = FALSE), "ggplot")
})
