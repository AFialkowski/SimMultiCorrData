library("SimMultiCorrData")
context("Lower skurtosis boundary")

skip_on_cran()
tol <- 1e-5
W <- calc_theory("Weibull", c(3, 5))
C <- calc_theory(fx = function(x) 0.7 * dchisq(x, 2) + 0.3 * dchisq(x, 32),
                 lower = 0, upper = Inf)
set.seed(104)
cstart1 <- runif(50, min = 0.5, max = 1.35)
cstart2 <- runif(50, min = -0.15, max = 0.05)
lstart1 <- runif(50, min = 0, max = 10)
xstart1 <- cbind(cstart1, cstart2, lstart1)
set.seed(104)
cstart1 <- runif(50, min = 0, max = 2)
cstart2 <- runif(50, min = -1, max = 1)
cstart3 <- runif(50, min = -1, max = 1)
cstart4 <- runif(50, min = -0.025, max = 0.025)
cstart5 <- runif(50, min = -0.025, max = 0.025)
lstart1 <- runif(50, min = -5, max = 5)
lstart2 <- runif(50, min = -0.025, max = 0.025)
lstart3 <- runif(50, min = -0.025, max = 0.025)
lstart4 <- runif(50, min = -0.025, max = 0.025)
xstart2 <- cbind(cstart1, cstart2, cstart3, cstart4,
                cstart5, lstart1, lstart2, lstart3, lstart4)

test_that("works for Fleishman method with a symmetric distribution", {
  expect_equal(all.equal(calc_lower_skurt(method = "Fleishman",
    skews = 0)$Min[1, "skurtosis"], -1.151323, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(calc_lower_skurt(method = "Fleishman",
    skews = 0, Skurt = seq(1.11, 2, 0.01))$Min[1, "skurtosis"], 0.008676821,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method with an asymmetric distribution", {
  expect_equal(all.equal(calc_lower_skurt(method = "Fleishman",
    skews = W[3])$Min[1, "skurtosis"], -1.102764, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(calc_lower_skurt(method = "Fleishman",
    skews = W[3], Skurt = seq(1.14, 2, 0.001))$Min[1, "skurtosis"],
    0.04423636, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method given xstart", {
  expect_equal(all.equal(calc_lower_skurt(method = "Fleishman", skews = W[3],
    xstart = xstart1, Skurt = seq(1.14, 2, 0.001))$Min[1, "skurtosis"],
    0.04423636, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method", {
  expect_equal(all.equal(calc_lower_skurt(method = "Polynomial", skews = W[3],
    fifths = W[5], sixths = W[6])$Min[1, "skurtosis"], -0.3234388,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(calc_lower_skurt(method = "Polynomial", skews = W[3],
    fifths = W[5], sixths = W[6] + 0.15,
    Skurt = seq(0.07, 2, 0.001))$Min[1, "skurtosis"], -0.271618,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method with a Six vector", {
  expect_error(calc_lower_skurt(method = "Polynomial", skews = C[3],
    fifths = C[5], sixths = C[6], Six = seq(0.01, 0.1, 0.01)))
})

test_that("works for Polynomial method given xstart", {
  expect_equal(all.equal(calc_lower_skurt(method = "Polynomial", skews = W[3],
    fifths = W[5], sixths = W[6] + 0.15, xstart = xstart2,
    Skurt = seq(0.07, 2, 0.001))$Min[1, "skurtosis"], -0.271618,
    tolerance = tol, check.attributes = FALSE), TRUE)
})
