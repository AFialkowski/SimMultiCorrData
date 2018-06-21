library("SimMultiCorrData")
context("Correlation boundaries by correlation method 1")

skip_on_cran()
options(scipen = 999)
tol <- 1e-5
L <- calc_theory("Logistic", c(0, 1))
Six <- list(seq(1.7, 1.8, 0.01))
marginal <- list(0.3)
lam <- 0.5
size <- 2
prob <- 0.75
mu <- size * (1 - prob)/prob
Rey <- matrix(0.4, 4, 4)
diag(Rey) <- 1

test_that("works for 0 continuous, 1 ordinal, 1 Poisson, 1 NB", {
  expect_equal(all.equal(valid_corr(k_cat = 1, k_pois = 1, k_nb = 1,
    marginal = marginal, lam = lam, size = size,
    prob = prob, rho = Rey[1:3, 1:3])$L_rho[1, 2], -0.7921534,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(valid_corr(k_cat = 1, k_pois = 1, k_nb = 1,
    marginal = marginal, lam = lam, size = size,
    mu = mu, rho = Rey[1:3, 1:3])$L_rho[1, 2], -0.7921534,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 1 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 1,
    k_nb = 1, method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal,
    lam = lam, size = size, prob = prob, rho = Rey)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 1,
    k_nb = 1, method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal,
    lam = lam, size = size, mu = mu, rho = Rey)$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 0 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 0, k_pois = 1,
    k_nb = 1, method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], lam = lam, size = size, prob = prob,
    rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 0, k_pois = 1,
    k_nb = 1, method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], lam = lam, size = size, mu = mu,
    rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 1 ordinal, 0 Poisson,
          1 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 0,
    k_nb = 1, method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal,
    size = size, prob = prob, rho = Rey[1:3, 1:3])$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 0,
    k_nb = 1, method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal,
    size = size, mu = mu, rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 1 ordinal, 1 Poisson,
          0 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal,
    lam = lam, rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 1 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 1,
    k_nb = 1, method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, lam = lam, size = size,
    prob = prob, rho = Rey)$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 1,
    k_nb = 1, method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, lam = lam, size = size,
    mu = mu, rho = Rey)$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 0 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 0, k_pois = 1,
    k_nb = 1, method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    lam = lam, size = size, prob = prob,
    rho = Rey[1:3, 1:3])$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 0, k_pois = 1,
    k_nb = 1, method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    lam = lam, size = size, mu = mu, rho = Rey[1:3, 1:3])$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 1 ordinal, 0 Poisson,
          1 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 0,
    k_nb = 1, method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, size = size,
    prob = prob, rho = Rey[1:3, 1:3])$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 0,
    k_nb = 1, method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six, marginal = marginal,
    size = size, mu = mu, rho = Rey[1:3, 1:3])$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 1 ordinal, 1 Poisson,
          0 NB", {
  expect_equal(all.equal(valid_corr(k_cont = 1, k_cat = 1, k_pois = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, lam = lam,
    rho = Rey[1:3, 1:3])$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

Rey2 <- matrix(0.4, 5, 5)
diag(Rey2) <- 1

test_that("works for Polynomial method: same continuous distribution", {
  expect_equal(all.equal(valid_corr(k_cont = 2, k_cat = 1, k_pois = 1,
    k_nb = 1, method = "Polynomial", means = rep(L[1], 2),
    vars = rep(L[2]^2, 2), skews = rep(L[3], 2), skurts = rep(L[4], 2),
    fifths = rep(L[5], 2), sixths = rep(L[6], 2), Six = list(1.75, 1.75),
    marginal = marginal, lam = lam, size = size,
    prob = prob, rho = Rey2)$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

