library("SimMultiCorrData")
context("Plots of theoretical pdf")

skip_on_cran()
W <- calc_theory("Weibull", c(3, 5))
Six <- seq(0.1, 1, 0.01)
WF <- nonnormvar1(method = "Fleishman", means = W[1], vars = W[2]^2,
                  skews = W[3], skurts = W[4])
WP <- nonnormvar1(method = "Polynomial", means = W[1], vars = W[2]^2,
                  skews = W[3], skurts = W[4], fifths = W[5],
                  sixths = W[6], Six = Six)
CF <- WF$constants
CP <- WP$constants

test_that("works for Fleishman method", {
  expect_is(plot_pdf_theory(CF, "Fleishman", mu = W[1], sigma = W[2],
    overlay = FALSE), "ggplot")
  expect_is(plot_pdf_theory(CF, "Fleishman", mu = W[1], sigma = W[2],
    overlay = TRUE, Dist = "Weibull", params = c(3, 5)), "ggplot")
  expect_is(plot_pdf_theory(CF, "Fleishman", mu = W[1], sigma = W[2],
    overlay = TRUE,
    fx = function(x) (3/5) * (x/5)^2 * exp(-(x/5)^3),
    lower = 0, upper = Inf), "ggplot")
})

test_that("works for Polynomial method", {
  expect_is(plot_pdf_theory(CP, "Polynomial", mu = W[1], sigma = W[2],
    overlay = FALSE), "ggplot")
  expect_is(plot_pdf_theory(CP, "Polynomial", mu = W[1], sigma = W[2],
    overlay = TRUE, Dist = "Weibull", params = c(3, 5)), "ggplot")
  expect_is(plot_pdf_theory(CP, "Polynomial", mu = W[1], sigma = W[2],
    overlay = TRUE,
    fx = function(x) (3/5) * (x/5)^2 * exp(-(x/5)^3),
    lower = 0, upper = Inf), "ggplot")
})

test_that("works for Fleishman method plus external data", {
  expect_is(plot_pdf_ext(CF, "Fleishman",
    ext_y = WF$continuous_variable[, 1]), "ggplot")
})

test_that("works for Polynomial method plus external data", {
  expect_is(plot_pdf_ext(CP, "Polynomial",
    ext_y = WP$continuous_variable[, 1]), "ggplot")
})
