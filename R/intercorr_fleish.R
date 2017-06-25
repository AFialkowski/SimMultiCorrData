#' @title Fleishman's Third-Order Polynomial Transformation Intermediate Correlation Equations
#'
#' @description This function contains Fleishman's third-order polynomial transformation intermediate correlation
#'     equations.  It is used in \code{\link[SimMultiCorrData]{findintercorr}} and \code{\link[SimMultiCorrData]{findintercorr2}}
#'     to find the intermediate correlation for standard normal random variables which are used in the Fleishman
#'     polynomial transformation.  It can be used to verify a set of constants and an intermediate correlation satisfy
#'     the equations for the desired post-transformation correlation.  It works for two or three variables.  Headrick &
#'     Sawilowsky (1999) recommend using the technique of Vale & Maurelli (1983) in the case of more than 3 variables, in which
#'     the intermediate correlations are found pairwise and then eigen value decomposition is used on the intermediate
#'     correlation matrix.  Note that there exist solutions that yield invalid power
#'     method pdfs (see \code{\link[SimMultiCorrData]{power_norm_corr}}, \code{\link[SimMultiCorrData]{pdf_check}}).
#'     This function would not ordinarily be called by the user.
#' @param r either a scalar, in which case it represents pairwise intermediate correlation between standard normal variables,
#'     or a vector of 3 values, in which case:
#'     \deqn{r[1]*r[2] = \rho_{z1,z2},\ r[1]*r[3] = \rho_{z1,z3},\ r[2]*r[3] = \rho_{z2,z3}}
#' @param c a matrix with either 2 or 3 rows, each a vector of constants c0, c1, c2, c3, like that returned by
#'     \code{\link[SimMultiCorrData]{find_constants}}
#' @param a a matrix of target correlations among continuous variables; if \code{nrow(a) = 1}, it represents a pairwise
#'     correlation; if \code{nrow(a) = 2 or 3}, it represents a correlation matrix between two or three variables
#' @export
#' @keywords intermediate, correlation, Fleishman
#' @seealso \code{\link[SimMultiCorrData]{fleish}}, \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{pdf_check}}, \code{\link[SimMultiCorrData]{find_constants}}
#' @return a list of length 1 for pairwise correlations or length 3 for three variables;
#'      if the inputs satisfy the equations, returns 0 for all list elements
#' @references Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#' Method. Psychometrika, 64, 25-35.
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#' Journal of Modern Applied Statistical Methods, 3, 65-71.
#'
#' Vale CD, Maurelli VA (1983). Simulating Multivariate Nonnormal Distributions. Psychometrika, 48, 465-471.
intercorr_fleish <- function(r, c, a) {
  if (nrow(a) == 1 | nrow(a) == 2) {
    f < -  numeric(1)
    a1 <- c[1, 1]
    b1 <- c[1, 2]
    d1 <- c[1, 4]
    a2 <- c[2, 1]
    b2 <- c[2, 2]
    d2 <- c[2, 4]
    if (nrow(a) == 1) rho <-  a[1, 1]
    if (nrow(a) == 2) rho <-  a[1, 2]
    f[1] <- r * (b1 * b2 + 3 * b2 * d1 + 3 * b1 * d2 + 9 * d1 * d2 +
                   2 * a1 * a2 * r + 6 * d1 * d2 * r^2) - rho
    return(f)
  }
  if (nrow(a) == 3) {
    f < -  numeric(3)
    a1 <- c[1, 1]
    b1 <- c[1, 2]
    d1 <- c[1, 4]
    a2 <- c[2, 1]
    b2 <- c[2, 2]
    d2 <- c[2, 4]
    a3 <- c[3, 1]
    b3 <- c[3, 2]
    d3 <- c[3, 4]
    rho1 <-  a[1, 2]
    rho2 <-  a[1, 3]
    rho3 <-  a[2, 3]
    f[1] <- r[1] * r[2] * (b1 * b2 + 3 * b2 * d1 + 3 * b1 * d2 +
                             9 * d1 * d2 + 2 * a1 * a2 * r[1] * r[2] +
                             6 * d1 * d2 * (r[1] * r[2])^2) - rho1
    f[2] <- r[1] * r[3] * (b1 * b3 + 3 * b3 * d1 + 3 * b1 * d3 +
                             9 * d1 * d3 + 2 * a1 * a3 * r[1] * r[3] +
                             6 * d1 * d3 * (r[1] * r[3])^2) - rho2
    f[3] <- r[2] * r[3] * (b2 * b3 + 3 * b3 * d2 + 3 * b2 * d3 +
                             9 * d2 * d3 + 2 * a2 * a3 * r[2] * r[3] +
                             6 * d2 * d3 * (r[2] * r[3])^2) - rho3
    return(f)
  }
}
