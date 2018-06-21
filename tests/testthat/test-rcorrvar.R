library("SimMultiCorrData")
context("Simulate using correlation method 1")

skip_on_cran()
options(scipen = 999)
tol <- 1e-5
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
  expect_equal(all.equal(rcorrvar(k_cat = 1, k_pois = 1, k_nb = 1,
    marginal = marginal, support = list(c(0, 1)), lam = lam, size = size,
    prob = prob, rho = Rey[1:3, 1:3])$maxerr, 0.06564162,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cat = 1, k_pois = 1, k_nb = 1,
    marginal = marginal, support = list(c(0, 1)), lam = lam, size = size,
    mu = mu, rho = Rey[1:3, 1:3])$maxerr, 0.06564162,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cat = 1, k_pois = 1, k_nb = 1,
    marginal = marginal, support = list(c(0, 1)), lam = lam, size = size,
    prob = prob, rho = Rey[1:3, 1:3], errorloop = TRUE)$maxerr, 0.0054317,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 1 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    lam = lam, size = size, prob = prob, rho = Rey)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    lam = lam, size = size, mu = mu, rho = Rey,
    cstart = list(cstartF))$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    lam = lam, size = size, prob = prob, rho = Rey,
    errorloop = TRUE)$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 0 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 0, k_pois = 1, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], lam = lam, size = size, prob = prob,
    rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 0, k_pois = 1, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], lam = lam, size = size, mu = mu,
    rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 0, k_pois = 1, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], lam = lam, size = size, prob = prob, rho = Rey[1:3, 1:3],
    errorloop = TRUE)$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 1 ordinal, 0 Poisson,
          1 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 0, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    size = size, prob = prob, rho = Rey[1:3, 1:3])$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 0, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    size = size, mu = mu, rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 0, k_nb = 1,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    size = size, prob = prob, rho = Rey[1:3, 1:3],
    errorloop = TRUE)$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 continuous, 1 ordinal, 1 Poisson,
          0 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 0,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    lam = lam, rho = Rey[1:3, 1:3])$constants[1, "c3"], 0.03605955,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 0,
    method = "Fleishman", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], marginal = marginal, support = list(c(0, 1)),
    lam = lam, rho = Rey[1:3, 1:3], errorloop = TRUE)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 1 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, support = list(c(0, 1)), lam = lam, size = size,
    prob = prob, rho = Rey)$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    cstart = list(cstartP), marginal = marginal, support = list(c(0, 1)),
    lam = lam, size = size, mu = mu, rho = Rey)$constants[1, "c5"],
    0.0000006125703, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, support = list(c(0, 1)), lam = lam, size = size,
    prob = prob, rho = Rey, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 0 ordinal, 1 Poisson,
          1 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 0, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    lam = lam, size = size, prob = prob,
    rho = Rey[1:3, 1:3])$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 0, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    lam = lam, size = size, mu = mu, rho = Rey[1:3, 1:3])$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 0, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    lam = lam, size = size, prob = prob, rho = Rey[1:3, 1:3],
    errorloop = TRUE)$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 1 ordinal, 0 Poisson,
          1 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 0, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, support = list(c(0, 1)), size = size,
    prob = prob, rho = Rey[1:3, 1:3])$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 0, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, support = list(c(0, 1)),
    size = size, mu = mu, rho = Rey[1:3, 1:3])$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 0, k_nb = 1,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, support = list(c(0, 1)), size = size,
    prob = prob, rho = Rey[1:3, 1:3], errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 continuous, 1 ordinal, 1 Poisson,
          0 NB", {
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 0,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, support = list(c(0, 1)), lam = lam,
    rho = Rey[1:3, 1:3])$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 1, k_cat = 1, k_pois = 1, k_nb = 0,
    method = "Polynomial", means = L[1], vars = L[2]^2, skews = L[3],
    skurts = L[4], fifths = L[5], sixths = L[6], Six = Six,
    marginal = marginal, support = list(c(0, 1)), lam = lam,
    rho = Rey[1:3, 1:3], errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

Rey2 <- matrix(0.4, 5, 5)
diag(Rey2) <- 1

test_that("works for Polynomial method: same continuous distribution", {
  expect_equal(all.equal(rcorrvar(k_cont = 2, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = rep(L[1], 2), vars = rep(L[2]^2, 2),
    skews = rep(L[3], 2), skurts = rep(L[4], 2), fifths = rep(L[5], 2),
    sixths = rep(L[6], 2), Six = list(1.75, 1.75),
    marginal = marginal, support = list(c(0, 1)), lam = lam, size = size,
    prob = prob, rho = Rey2)$constants[1, "c5"], 0.0000006124845,
    tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(rcorrvar(k_cont = 2, k_cat = 1, k_pois = 1, k_nb = 1,
    method = "Polynomial", means = rep(L[1], 2), vars = rep(L[2]^2, 2),
    skews = rep(L[3], 2), skurts = rep(L[4], 2), fifths = rep(L[5], 2),
    sixths = rep(L[6], 2), Six = list(1.75, 1.75),
    marginal = marginal, support = list(c(0, 1)), lam = lam, size = size,
    prob = prob, rho = Rey2, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

