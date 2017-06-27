#' @title Fleishman's Third-Order Polynomial Transformation Equations
#'
#' @description This function contains Fleishman's third-order polynomial transformation equations (\doi{10.1007/BF02293811}).  It is used in
#'     \code{\link[SimMultiCorrData]{find_constants}} to find the constants c1, c2, and c3 (c0 = -c2) that satisfy the
#'     equations given skewness and standardized kurtosis values.  It can be used to verify a set of constants satisfy
#'     the equations.  Note that there exist solutions that yield invalid power method pdfs (see
#'     \code{\link[SimMultiCorrData]{power_norm_corr}}, \code{\link[SimMultiCorrData]{pdf_check}}).
#'     This function would not ordinarily be called by the user.
#' @param c a vector of constants c1, c2, c3; note that \code{\link[SimMultiCorrData]{find_constants}} returns c0, c1, c2, c3
#' @param a a vector c(skewness, standardized kurtosis)
#' @export
#' @keywords constants, Fleishman
#' @seealso \code{\link[SimMultiCorrData]{poly}}, \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{pdf_check}}, \code{\link[SimMultiCorrData]{find_constants}}
#' @return a list of length 3; if the constants satisfy the equations, returns 0 for all list elements
#' @references
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532. \doi{10.1007/BF02293811}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35. \doi{10.1007/BF02294317}.
#'
#' @examples
#' # Laplace Distribution
#' fleish(c = c(0.782356, 0, 0.067905), a = c(0, 3))
fleish <- function(c, a) {
  f <- numeric(3)
  g1 <- a[1]
  g2 <- a[2]
  f[1] <- c[1]^2 + 6 * c[1] * c[3] + 2 * c[2]^2 + 15 * c[3]^2 - 1
  f[2] <- 2 * c[2] * (c[1]^2 + 24 * c[1] * c[3] + 105 * c[3]^2 + 2) - g1
  f[3] <- 24 * (c[1] * c[3] + c[2]^2 * (1 + c[1]^2 + 28 * c[1] * c[3]) +
                  c[3]^2 * (12 + 48 * c[1] * c[3] + 141 * c[2]^2 +
                              225 * c[3]^2)) - g2
  return(f)
}
