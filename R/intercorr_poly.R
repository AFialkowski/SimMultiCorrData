#' @title Headrick's Fifth-Order Polynomial Transformation Intermediate Correlation Equations
#'
#' @description This function contains Headrick's fifth-order polynomial transformation intermediate correlation
#'     equations (2002, \doi{10.1016/S0167-9473(02)00072-5}).  It is used in \code{\link[SimMultiCorrData]{findintercorr}} and
#'     \code{\link[SimMultiCorrData]{findintercorr2}}
#'     to find the intermediate correlation for standard normal random variables which are used in the Headrick
#'     polynomial transformation.  It can be used to verify a set of constants and an intermediate correlation satisfy
#'     the equations for the desired post-transformation correlation.  It works for two, three, or four variables.  Headrick recommended
#'     using the technique of Vale & Maurelli (1983,
#'     \doi{10.1007/BF02293687}), in the case of more than 4 variables, in which
#'     the intermediate correlations are found pairwise and then eigen value decomposition is used on the correlation matrix.
#'     Note that there exist solutions that yield invalid power
#'     method pdfs (see \code{\link[SimMultiCorrData]{power_norm_corr}}, \code{\link[SimMultiCorrData]{pdf_check}}).
#'     This function would not ordinarily be called by the user.
#' @param r either a scalar, in which case it represents pairwise intermediate correlation between standard normal variables,
#'     or a vector of 3 values, in which case:
#'     \deqn{r[1]*r[2] = \rho_{z1,z2},\ r[1]*r[3] = \rho_{z1,z3},\ r[2]*r[3] = \rho_{z2,z3}}
#'     or a vector of 4 values, in which case:
#'     \deqn{r0 = r[5]*r[6],\ r0*r[1]*r[2] = \rho_{z1,z2},\ r0*r[1]*r[3] = \rho_{z1,z3}}
#'     \deqn{r0*r[2]*r[3] = \rho_{z2,z3},\ r0*r[1]*r[4] = \rho_{z1,z4},\ r0*r[2]*r[4] = \rho_{z2,z4},}
#'     \deqn{r0*r[3]*r[4] = \rho_{z3,z4}}
#' @param c a matrix with either 2, 3, or 4 rows, each a vector of constants c0, c1, c2, c3, like that returned by
#'     \code{\link[SimMultiCorrData]{find_constants}}
#' @param a a matrix of target correlations among continuous variables; if \code{nrow(a) = 1}, it represents a pairwise
#'     correlation; if \code{nrow(a) = 2, 3, or 4}, it represents a correlation matrix between two, three, or four variables
#' @export
#' @keywords intermediate, correlation, Headrick
#' @seealso \code{\link[SimMultiCorrData]{poly}}, \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{pdf_check}}, \code{\link[SimMultiCorrData]{find_constants}}
#' @return a list of length 1 for pairwise correlations, length 3 for three variables, or length 6 for four variables;
#'      if the inputs satisfy the equations, returns 0 for all list elements
#' @references Please see references for \code{\link[SimMultiCorrData]{findintercorr_cont}}.
#'
intercorr_poly <- function(r, c, a) {
  if (nrow(a) == 1 | nrow(a) == 2) {
    f <- numeric(1)
    c0a <- c[1, 1]
    c1a <- c[1, 2]
    c2a <- c[1, 3]
    c3a <- c[1, 4]
    c4a <- c[1, 5]
    c5a <- c[1, 6]
    c0b <- c[2, 1]
    c1b <- c[2, 2]
    c2b <- c[2, 3]
    c3b <- c[2, 4]
    c4b <- c[2, 5]
    c5b <- c[2, 6]
    if (nrow(a) == 1) rho <- a[1, 1]
    if (nrow(a) == 2) rho <- a[1, 2]
    f[1] <- (3 * c4a * c0b + 3 * c4a * c2b + 9 * c4a * c4b + c0a *
               (c0b + c2b + 3 * c4b) + c1a * c1b * r + 3 * c3a * c1b * r +
               15 * c5a * c1b * r + 3 * c1a * c3b * r + 9 * c3a * c3b * r +
               45 * c5a * c3b * r + 15 * c1a * c5b * r + 45 * c3a * c5b * r +
               225 * c5a * c5b * r + 12 * c4a * c2b * r^2 +
               72 * c4a * c4b * r^2 + 6 * c3a * c3b * r^3 +
               60 * c5a * c3b * r^3 + 60 * c3a * c5b * r^3 +
               600 * c5a * c5b * r^3 + 24 * c4a * c4b * r^4 +
               120 * c5a * c5b * r^5 + c2a * (c0b + c2b + 3 * c4b +
               2 * c2b * r^2 + 12 * c4b * r^2)) - rho
    return(f)
  }
  if (nrow(a) == 3) {
    f <- numeric(3)
    c0a <- c[1, 1]
    c1a <- c[1, 2]
    c2a <- c[1, 3]
    c3a <- c[1, 4]
    c4a <- c[1, 5]
    c5a <- c[1, 6]
    c0b <- c[2, 1]
    c1b <- c[2, 2]
    c2b <- c[2, 3]
    c3b <- c[2, 4]
    c4b <- c[2, 5]
    c5b <- c[2, 6]
    c0c <- c[3, 1]
    c1c <- c[3, 2]
    c2c <- c[3, 3]
    c3c <- c[3, 4]
    c4c <- c[3, 5]
    c5c <- c[3, 6]
    rho1 <- a[1, 2]
    rho2 <- a[1, 3]
    rho3 <- a[2, 3]
    f[1] <- (3 * c4a * c0b + 3 * c4a * c2b + 9 * c4a * c4b +
               c0a * (c0b + c2b + 3 * c4b) + c1a * c1b * r[1] * r[2] +
               3 * c3a * c1b * r[1] * r[2] + 15 * c5a * c1b * r[1] * r[2] +
               3 * c1a * c3b * r[1] * r[2] + 9 * c3a * c3b * r[1] * r[2] +
               45 * c5a * c3b * r[1] * r[2] + 15 * c1a * c5b * r[1] * r[2] +
               45 * c3a * c5b * r[1] * r[2] + 225 * c5a * c5b * r[1] * r[2] +
               12 * c4a * c2b * (r[1] * r[2])^2 + 72 * c4a * c4b *
               (r[1] * r[2])^2 + 6 * c3a * c3b * (r[1] * r[2])^3 +
               60 * c5a * c3b * (r[1] * r[2])^3 + 60 * c3a * c5b *
               (r[1] * r[2])^3 + 600 * c5a * c5b * (r[1] * r[2])^3 +
               24 * c4a * c4b * (r[1] * r[2])^4 + 120 * c5a * c5b *
               (r[1] * r[2])^5 + c2a * (c0b + c2b + 3 * c4b + 2 * c2b *
              (r[1] * r[2])^2 + 12 * c4b * (r[1] * r[2])^2)) - rho1
    f[2] <- (3 * c4a * c0c + 3 * c4a * c2c + 9 * c4a * c4c + c0a *
               (c0c + c2c + 3 * c4c) + c1a * c1c * r[1] * r[3] +
               3 * c3a * c1c * r[1] * r[3] + 15 * c5a * c1c * r[1] * r[3] +
               3 * c1a * c3c * r[1] * r[3] + 9 * c3a * c3c * r[1] * r[3] +
               45 * c5a * c3c * r[1] * r[3] + 15 * c1a * c5c * r[1] * r[3] +
               45 * c3a * c5c * r[1] * r[3] + 225 * c5a * c5c * r[1] * r[3] +
               12 * c4a * c2c * (r[1] * r[3])^2 + 72 * c4a * c4c *
               (r[1] * r[3])^2 + 6 * c3a * c3c * (r[1] * r[3])^3 +
               60 * c5a * c3c * (r[1] * r[3])^3 + 60 * c3a * c5c *
               (r[1] * r[3])^3 + 600 * c5a * c5c * (r[1] * r[3])^3 +
               24 * c4a * c4c * (r[1] * r[3])^4 + 120 * c5a * c5c *
               (r[1] * r[3])^5 + c2a * (c0c + c2c + 3 * c4c + 2 * c2c *
              (r[1] * r[3])^2 + 12 * c4c * (r[1] * r[3])^2)) - rho2
    f[3] <- (3 * c4b * c0c + 3 * c4b * c2c + 9 * c4b * c4c + c0b *
               (c0c + c2c + 3 * c4c) + c1b * c1c * r[2] * r[3] +
               3 * c3b * c1c * r[2] * r[3] + 15 * c5b * c1c * r[2] * r[3] +
               3 * c1b * c3c * r[2] * r[3] + 9 * c3b * c3c * r[2] * r[3] +
               45 * c5b * c3c * r[2] * r[3] + 15 * c1b * c5c * r[2] * r[3] +
               45 * c3b * c5c * r[2] * r[3] + 225 * c5b * c5c * r[2] * r[3] +
               12 * c4b * c2c * (r[2] * r[3])^2 + 72 * c4b * c4c *
               (r[2] * r[3])^2 + 6 * c3b * c3c * (r[2] * r[3])^3 +
               60 * c5b * c3c * (r[2] * r[3])^3 + 60 * c3b * c5c *
               (r[2] * r[3])^3 + 600 * c5b * c5c * (r[2] * r[3])^3 +
               24 * c4b * c4c * (r[2] * r[3])^4 + 120 * c5b * c5c *
               (r[2] * r[3])^5 + c2b * (c0c + c2c + 3 * c4c + 2 * c2c *
               (r[2] * r[3])^2 + 12 * c4c * (r[2] * r[3])^2)) - rho3
    return(f)
  }
  if (nrow(a) == 4) {
    f <- numeric(6)
    c0a <- c[1, 1]
    c1a <- c[1, 2]
    c2a <- c[1, 3]
    c3a <- c[1, 4]
    c4a <- c[1, 5]
    c5a <- c[1, 6]
    c0b <- c[2, 1]
    c1b <- c[2, 2]
    c2b <- c[2, 3]
    c3b <- c[2, 4]
    c4b <- c[2, 5]
    c5b <- c[2, 6]
    c0c <- c[3, 1]
    c1c <- c[3, 2]
    c2c <- c[3, 3]
    c3c <- c[3, 4]
    c4c <- c[3, 5]
    c5c <- c[3, 6]
    c0d <- c[4, 1]
    c1d <- c[4, 2]
    c2d <- c[4, 3]
    c3d <- c[4, 4]
    c4d <- c[4, 5]
    c5d <- c[4, 6]
    rho1 <- a[1, 2]
    rho2 <- a[1, 3]
    rho3 <- a[2, 3]
    rho4 <- a[1, 4]
    rho5 <- a[2, 4]
    rho6 <- a[3, 4]
    f[1] <- (3 * c4a * c0b + 3 * c4a * c2b + 9 * c4a * c4b + c0a *
               (c0b + c2b + 3 * c4b) + c1a * c1b * r[5] * r[6] * r[1] * r[2] +
               3 * c3a * c1b * r[5] * r[6] * r[1] * r[2] + 15 * c5a * c1b *
               r[5] * r[6] * r[1] * r[2] + 3 * c1a * c3b * r[5] * r[6] *
               r[1] * r[2] + 9 * c3a * c3b * r[5] * r[6] * r[1] * r[2] +
               45 * c5a * c3b * r[5] * r[6] * r[1] * r[2] + 15 * c1a * c5b *
               r[5] * r[6] * r[1] * r[2] + 45 * c3a * c5b * r[5] * r[6] *
               r[1] * r[2] + 225 * c5a * c5b * r[5] * r[6] * r[1] * r[2] +
               12 * c4a * c2b * (r[5] * r[6] * r[1] * r[2])^2 + 72 * c4a *
               c4b * (r[5] * r[6] * r[1] * r[2])^2 + 6 * c3a * c3b * (r[5] *
               r[6] * r[1] * r[2])^3 + 60 * c5a * c3b * (r[5] * r[6] * r[1] *
               r[2])^3 + 60 * c3a * c5b * (r[5] * r[6] * r[1] * r[2])^3 +
               600 * c5a * c5b * (r[5] * r[6] * r[1] * r[2])^3 + 24 * c4a *
               c4b * (r[5] * r[6] * r[1] * r[2])^4 + 120 * c5a * c5b *
               (r[5] * r[6] * r[1] * r[2])^5 + c2a * (c0b + c2b + 3 * c4b +
               2 * c2b * (r[5] * r[6] * r[1] * r[2])^2 + 12 * c4b * (r[5] *
               r[6] * r[1] * r[2])^2)) - rho1
    f[2] <- (3 * c4a * c0c + 3 * c4a * c2c + 9 * c4a * c4c + c0a *
               (c0c + c2c + 3 * c4c) + c1a * c1c * r[5] * r[6] * r[1] * r[3] +
               3 * c3a * c1c * r[5] * r[6] * r[1] * r[3] + 15 * c5a * c1c *
               r[5] * r[6] * r[1] * r[3] + 3 * c1a * c3c * r[5] * r[6] *
               r[1] * r[3] + 9 * c3a * c3c * r[5] * r[6] * r[1] * r[3] +
               45 * c5a * c3c * r[5] * r[6] * r[1] * r[3] + 15 * c1a * c5c *
               r[5] * r[6] * r[1] * r[3] + 45 * c3a * c5c * r[5] * r[6] *
               r[1] * r[3] + 225 * c5a * c5c * r[5] * r[6] * r[1] * r[3] +
               12 * c4a * c2c * (r[5] * r[6] * r[1] * r[3])^2 + 72 * c4a *
               c4c * (r[5] * r[6] * r[1] * r[3])^2 + 6 * c3a * c3c * (r[5] *
               r[6] * r[1] * r[3])^3 + 60 * c5a * c3c * (r[5] * r[6] * r[1] *
               r[3])^3 + 60 * c3a * c5c * (r[5] * r[6] * r[1] * r[3])^3 +
               600 * c5a * c5c * (r[5] * r[6] * r[1] * r[3])^3 + 24 * c4a *
               c4c * (r[5] * r[6] * r[1] * r[3])^4 + 120 * c5a * c5c * (r[5] *
               r[6] * r[1] * r[3])^5 + c2a * (c0c + c2c + 3 * c4c + 2 * c2c *
               (r[5] * r[6] * r[1] * r[3])^2 + 12 * c4c * (r[5] * r[6] *
               r[1] * r[3])^2)) - rho2
    f[3] <- (3 * c4b * c0c + 3 * c4b * c2c + 9 * c4b * c4c + c0b *
               (c0c + c2c + 3 * c4c) + c1b * c1c * r[5] * r[6] * r[2] * r[3] +
               3 * c3b * c1c * r[5] * r[6] * r[2] * r[3] + 15 * c5b * c1c *
               r[5] * r[6] * r[2] * r[3] + 3 * c1b * c3c * r[5] * r[6] *
               r[2] * r[3] + 9 * c3b * c3c * r[5] * r[6] * r[2] * r[3] +
               45 * c5b * c3c * r[5] * r[6] * r[2] * r[3] + 15 * c1b * c5c *
               r[5] * r[6] * r[2] * r[3] + 45 * c3b * c5c * r[5] * r[6] *
               r[2] * r[3] + 225 * c5b * c5c * r[5] * r[6] * r[2] * r[3] +
               12 * c4b * c2c * (r[5] * r[6] * r[2] * r[3])^2 + 72 * c4b *
               c4c * (r[5] * r[6] * r[2] * r[3])^2 + 6 * c3b * c3c * (r[5] *
               r[6] * r[2] * r[3])^3 + 60 * c5b * c3c * (r[5] * r[6] *
               r[2] * r[3])^3 + 60 * c3b * c5c * (r[5] * r[6] * r[2] *
               r[3])^3 + 600 * c5b * c5c * (r[5] * r[6] * r[2] * r[3])^3 +
               24 * c4b * c4c * (r[5] * r[6] * r[2] * r[3])^4 + 120 * c5b *
               c5c * (r[5] * r[6] * r[2] * r[3])^5 + c2b * (c0c + c2c + 3 *
               c4c + 2 * c2c * (r[5] * r[6] * r[2] * r[3])^2 + 12 * c4c *
              (r[5] * r[6] * r[2] * r[3])^2)) - rho3
    f[4] <- (3 * c4a * c0d + 3 * c4a * c2d + 9 * c4a * c4d + c0a *
               (c0d + c2d + 3 * c4d) + c1a * c1d * r[5] * r[6] * r[1] * r[4] +
               3 * c3a * c1d * r[5] * r[6] * r[1] * r[4] + 15 * c5a * c1d *
               r[5] * r[6] * r[1] * r[4] + 3 * c1a * c3d * r[5] * r[6] *
               r[1] * r[4] + 9 * c3a * c3d * r[5] * r[6] * r[1] * r[4] +
               45 * c5a * c3d * r[5] * r[6] * r[1] * r[4] + 15 * c1a *
               c5d * r[5] * r[6] * r[1] * r[4] + 45 * c3a * c5d * r[5] *
               r[6] * r[1] * r[4] + 225 * c5a * c5d * r[5] * r[6] * r[1] *
               r[4] + 12 * c4a * c2d * (r[5] * r[6] * r[1] * r[4])^2 + 72 *
               c4a * c4d * (r[5] * r[6] * r[1] * r[4])^2 + 6 * c3a * c3d *
               (r[5] * r[6] * r[1] * r[4])^3 + 60 * c5a * c3d * (r[5] *
               r[6] * r[1] * r[4])^3 + 60 * c3a * c5d * (r[5] * r[6] *
               r[1] * r[4])^3 + 600 * c5a * c5d * (r[5] * r[6] * r[1] *
               r[4])^3 + 24 * c4a * c4d * (r[5] * r[6] * r[1] * r[4])^4 +
               120 * c5a * c5d * (r[5] * r[6] * r[1] * r[4])^5 + c2a *
               (c0d + c2d + 3 * c4d + 2 * c2d * (r[5] * r[6] * r[1] *
              r[4])^2 + 12 * c4d * (r[5] * r[6] * r[1] * r[4])^2)) - rho4
    f[5] <- (3 * c4b * c0d + 3 * c4b * c2d + 9 * c4b * c4d + c0b *
               (c0d + c2d + 3 * c4d) + c1b * c1d * r[5] * r[6] * r[2] * r[4] +
               3 * c3b * c1d * r[5] * r[6] * r[2] * r[4] + 15 * c5b * c1d *
               r[5] * r[6] * r[2] * r[4] + 3 * c1b * c3d * r[5] * r[6] *
               r[2] * r[4] + 9 * c3b * c3d * r[5] * r[6] * r[2] * r[4] +
               45 * c5b * c3d * r[5] * r[6] * r[2] * r[4] + 15 * c1b * c5d *
               r[5] * r[6] * r[2] * r[4] + 45 * c3b * c5d * r[5] * r[6] *
               r[2] * r[4] + 225 * c5b * c5d * r[5] * r[6] * r[2] * r[4] +
               12 * c4b * c2d * (r[5] * r[6] * r[2] * r[4])^2 + 72 * c4b *
               c4d * (r[5] * r[6] * r[2] * r[4])^2 + 6 * c3b * c3d *
               (r[5] * r[6] * r[2] * r[4])^3 + 60 * c5b * c3d * (r[5] *
               r[6] * r[2] * r[4])^3 + 60 * c3b * c5d * (r[5] * r[6] *
               r[2] * r[4])^3 + 600 * c5b * c5d * (r[5] * r[6] * r[2] *
               r[4])^3 + 24 * c4b * c4d * (r[5] * r[6] * r[2] * r[4])^4 +
               120 * c5b * c5d * (r[5] * r[6] * r[2] * r[4])^5 + c2b *
               (c0d + c2d + 3 * c4d + 2 * c2d * (r[5] * r[6] * r[2] *
              r[4])^2 + 12 * c4d * (r[5] * r[6] * r[2] * r[4])^2)) - rho5
    f[6] <- (3 * c4c * c0d + 3 * c4c * c2d + 9 * c4c * c4d + c0c *
               (c0d + c2d + 3 * c4d) + c1c * c1d * r[5] * r[6] * r[3] * r[4] +
               3 * c3c * c1d * r[5] * r[6] * r[3] * r[4] + 15 * c5c * c1d *
               r[5] * r[6] * r[3] * r[4] + 3 * c1c * c3d * r[5] * r[6] *
               r[3] * r[4] + 9 * c3c * c3d * r[5] * r[6] * r[3] * r[4] +
               45 * c5c * c3d * r[5] * r[6] * r[3] * r[4] + 15 * c1c * c5d *
               r[5] * r[6] * r[3] * r[4] + 45 * c3c * c5d * r[5] * r[6] *
               r[3] * r[4] + 225 * c5c * c5d * r[5] * r[6] * r[3] * r[4] +
               12 * c4c * c2d * (r[5] * r[6] * r[3] * r[4])^2 + 72 * c4c *
               c4d * (r[5] * r[6] * r[3] * r[4])^2 + 6 * c3c * c3d *
               (r[5] * r[6] * r[3] * r[4])^3 + 60 * c5c * c3d * (r[5] *
               r[6] * r[3] * r[4])^3 + 60 * c3c * c5d * (r[5] * r[6] * r[3] *
               r[4])^3 + 600 * c5c * c5d * (r[5] * r[6] * r[3] * r[4])^3 +
               24 * c4c * c4d * (r[5] * r[6] * r[3] * r[4])^4 + 120 * c5c *
               c5d * (r[5] * r[6] * r[3] * r[4])^5 + c2c * (c0d + c2d +
               3 * c4d + 2 * c2d * (r[5] * r[6] * r[3] * r[4])^2 +
              12 * c4d * (r[5] * r[6] * r[3] * r[4])^2)) - rho6
    return(f)
  }
}
