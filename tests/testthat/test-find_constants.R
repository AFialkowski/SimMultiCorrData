library("SimMultiCorrData")
context("Power method constants")

skip_on_cran()
options(scipen = 999)
tol <- 1e-5
L <- calc_theory("Logistic", c(0, 1))
W <- calc_theory("Weibull", c(3, 5))
set.seed(1234)
n <- 25
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -0.5, max = 0.5)
cstartF <- cbind(cstart1, cstart2, cstart3)
set.seed(1234)
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -1, max = 1)
cstart4 <- runif(n, min = -0.025, max = 0.025)
cstart5 <- runif(n, min = -0.025, max = 0.025)
cstartP <- cbind(cstart1, cstart2, cstart3, cstart4, cstart5)

test_that("works for Fleishman method with a symmetric distribution", {
  expect_equal(all.equal(find_constants(method = "Fleishman", skews = L[3],
    skurts = L[4])$constants["c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method with an asymmetric distribution", {
  expect_equal(all.equal(find_constants(method = "Fleishman", skews = W[3],
    skurts = W[4])$constants["c3"], -0.01467219, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method given cstart", {
  expect_equal(all.equal(find_constants(method = "Fleishman", skews = W[3],
    skurts = W[4], cstart = cstartF)$constants["c3"], -0.01467219,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method", {
  expect_equal(all.equal(find_constants(method = "Polynomial", skews = W[3],
    skurts = W[4], fifths = W[5], sixths = W[6])$constants["c5"],
    -0.002740699, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method with a Six vector", {
  expect_equal(all.equal(find_constants(method = "Polynomial", skews = W[3],
    skurts = W[4], fifths = W[5], sixths = W[6],
    Six = seq(0.01, 2, 0.01))$constants["c5"], 0.0004771859, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method given cstart", {
  expect_equal(all.equal(find_constants(method = "Polynomial", skews = W[3],
    skurts = W[4], fifths = W[5], sixths = W[6], cstart = cstartP,
    Six = seq(0.01, 2, 0.01))$constants["c5"], 0.0004771859, tolerance = tol,
    check.attributes = FALSE), TRUE)
})
