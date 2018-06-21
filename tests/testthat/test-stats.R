library("SimMultiCorrData")
context("Summary statistics")

skip_on_cran()
tol <- 1e-5
W <- calc_theory("Logistic", c(0, 1))
Six <- 1.75
CF <- find_constants("Fleishman", W[3], W[4])$constants
CP <- find_constants("Polynomial", W[3], W[4], W[5], W[6], Six = Six)$constants

test_that("pdf_check works for Fleishman method", {
  expect_equal(pdf_check(CF, "Fleishman")$valid.pdf, TRUE)
})

test_that("pdf_check works for Polynomial method", {
  expect_equal(pdf_check(CP, "Polynomial")$valid.pdf, TRUE)
})

test_that("corr works for Fleishman method", {
  expect_equal(all.equal(power_norm_corr(CF, "Fleishman"), 0.9960915,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("corr works for Polynomial method", {
  expect_equal(all.equal(power_norm_corr(CP, "Polynomial"), 0.9960917,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("stats works for Fleishman method", {
  expect_equal(all.equal(stats_pdf(CF, "Fleishman", mu = W[1],
    sigma = W[2])["mode"], 0, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("stats works for Polynomial method", {
  expect_equal(all.equal(stats_pdf(CP, "Polynomial", mu = W[1],
    sigma = W[2])["mode"], 0, tolerance = tol, check.attributes = FALSE), TRUE)
})
