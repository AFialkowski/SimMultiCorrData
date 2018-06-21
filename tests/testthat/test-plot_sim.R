library("SimMultiCorrData")
context("Plots of simulated data")

skip_on_cran()
W <- calc_theory("Weibull", c(3, 5))
Six <- seq(0.1, 1, 0.01)
WF <- nonnormvar1(method = "Fleishman", means = W[1], vars = W[2]^2,
                  skews = W[3], skurts = W[4])
WP <- nonnormvar1(method = "Polynomial", means = W[1], vars = W[2]^2,
                  skews = W[3], skurts = W[4], fifths = W[5],
                  sixths = W[6], Six = Six)

marginal <- list(0.3)
lam <- 5
size <- 2
prob <- 0.75
Rey <- matrix(0.4, 3, 3)
diag(Rey) <- 1
Sim1 <- rcorrvar(k_cat = 1, k_pois = 1, k_nb = 1, marginal = marginal,
  support = list(c(0, 1)), lam = lam, size = size, prob = prob, rho = Rey)

test_that("works for continuous variable", {
  expect_is(plot_sim_theory(WP$continuous_variable[, 1],
    overlay = FALSE), "ggplot")
  expect_is(plot_sim_theory(WP$continuous_variable[, 1],
    overlay = TRUE, Dist = "Weibull", params = c(3, 5)), "ggplot")
  expect_is(plot_sim_theory(WP$continuous_variable[, 1],
    overlay = TRUE,
    fx = function(x) (3/5) * (x/5)^2 * exp(-(x/5)^3),
    lower = 0.00001, upper = 1000), "ggplot")
})

test_that("works for Poisson variable", {
  expect_is(plot_sim_theory(Sim1$Poisson_variables[, 1],
    overlay = FALSE), "ggplot")
  expect_is(plot_sim_theory(Sim1$Poisson_variables[, 1], cont_var = FALSE,
    Dist = "Poisson", params = lam), "ggplot")
})

test_that("works for NB variable", {
  expect_is(plot_sim_theory(Sim1$Neg_Bin_variables[, 1],
    overlay = FALSE), "ggplot")
  expect_is(plot_sim_theory(Sim1$Neg_Bin_variables[, 1], cont_var = FALSE,
    Dist = "Negative_Binomial", params = c(size, prob)), "ggplot")
})

test_that("works for continuous variable plus external data", {
  expect_is(plot_sim_ext(WP$continuous_variable[, 1],
    ext_y = rweibull(10000, 3, 5)), "ggplot")
})

test_that("works for Poisson variable plus external data", {
  expect_is(plot_sim_ext(Sim1$Poisson_variables[, 1],
    ext_y = rpois(10000, lam)), "ggplot")
})

test_that("works for NB variable plus external data", {
  expect_is(plot_sim_ext(Sim1$Neg_Bin_variables[, 1],
    ext_y = rnbinom(10000, size, prob)), "ggplot")
})
