#' @title Fleishman Transformation Lagrangean Constraints for Lower Boundary of Standardized Kurtosis in Asymmetric Distributions
#'
#' @description This function gives the first-order conditions of the Fleishman Transformation Lagrangean expression
#'     \eqn{F(c1, c3, \lambda) = f(c1, c3) + \lambda * [\gamma_{1} - g(c1, c3)]} used to find the lower kurtosis boundary for a given non-zero skewness
#'     in \code{\link[SimMultiCorrData]{calc_lower_skurt}} (see Headrick & Sawilowsky, 2002).  Here, \eqn{f(c1, c3)} is the equation for
#'     standardized kurtosis expressed in terms of c1 and c3 only,
#'     \eqn{\lambda} is the Lagrangean multiplier, \eqn{\gamma_{1}} is skewness, and \eqn{g(c1, c3)} is the equation for skewness expressed
#'     in terms of c1 and c3 only.  It should be noted that these equations are for \eqn{\gamma_{1} > 0}.  Negative skew values are handled within
#'     \code{\link[SimMultiCorrData]{calc_lower_skurt}}.  Headrick & Sawilowsky (2002) gave equations for the first-order derivatives \eqn{dF/dc1}
#'     and \eqn{dF/dc3}.  These were verified and \eqn{dF/d\lambda} was calculated using \code{\link[stats]{D}}.  The second-order conditions to
#'     verify that the kurtosis is a global minimum are in \code{\link[SimMultiCorrData]{fleish_Hessian}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param c a vector of constants c1, c3, lambda
#' @param a skew value
#' @export
#' @keywords kurtosis, boundary, Fleishman
#' @seealso \code{\link[SimMultiCorrData]{fleish_Hessian}}, \code{\link[SimMultiCorrData]{calc_lower_skurt}}
#' @return A list with components:
#' @return \eqn{dF(c1, c3, \lambda)/d\lambda = \gamma_{1} - g(c1, c3)}
#' @return \eqn{dF(c1, c3, \lambda)/dc1 = df(c1, c3)/dc1  - \lambda * dg(c1, c3)/dc1}
#' @return \eqn{dF(c1, c3, \lambda)/dc3 = df(c1, c3)/dc3  - \lambda * dg(c1, c3)/dc3}
#' @return If the suppled values for c and skew satisfy the Lagrangean expression, it will return 0 for each component.
#' @references Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532.
#'
#' Headrick TC, Sawilowsky SS (2002). Weighted Simplex Procedures for Determining Boundary Points and Constants for the
#'     Univariate and Multivariate Power Methods. Journal of Educational and Behavioral Statistics, 25, 417-436.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35.
#'
fleish_skurt_check <- function(c, a) {
  F <- numeric(3)
  g1 <- a[1]
  F[1] <- g1 - ((8 - 6 * c[1]^4 - 2 * c[1]^6 + 144 * c[1] * c[2] -
                   144 * c[1]^3 * c[2] - 108 * c[1]^5 * c[2] +
                   720 * c[2]^2 - 540 * c[1]^2 * c[2]^2 -
                   2178 * c[1]^4 * c[2]^2 + 2160 * c[1] * c[2]^3 -
                   20952 * c[1]^3 * c[2]^3 + 9450 * c[2]^4 -
                   106110 * c[1]^2 * c[2]^4 - 283500 * c[1] * c[2]^5 -
                   330750 * c[2]^6)^(1/2))
  F[2] <- (12 * (24 * c[2] - 4 * c[1]^3 - 34 * (3 * c[1]^2) * c[2] -
                   324 * (2 * c[1]) * c[2]^2 - 1170 * c[2]^3)) -
    c[3] * ((8 - 6 * c[1]^4 - 2 * c[1]^6 + 144 * c[1] * c[2] -
               144 * c[1]^3 * c[2] - 108 * c[1]^5 * c[2] +
               720 * c[2]^2 - 540 * c[1]^2 * c[2]^2 -
               2178 * c[1]^4 * c[2]^2 + 2160 * c[1] * c[2]^3 -
               20952 * c[1]^3 * c[2]^3 + 9450 * c[2]^4 -
               106110 * c[1]^2 * c[2]^4 - 283500 * c[1] * c[2]^5 -
               330750 * c[2]^6)^(-1/2) *
              ((1/2) * (144 * c[2] - (6 * (4 * c[1]^3) +
                                        2 * (6 * c[1]^5)) -
                          144 * (3 * c[1]^2) * c[2] -
                          108 * (5 * c[1]^4) * c[2] -
                          540 * (2 * c[1]) * c[2]^2 -
                          2178 * (4 * c[1]^3) * c[2]^2 +
                          2160 * c[2]^3 -
                          20952 * (3 * c[1]^2) * c[2]^3 -
                          106110 * (2 * c[1]) * c[2]^4 -
                          283500 * c[2]^5)))
  F[3] <- (12 * (24 * c[1] - 34 * c[1]^3 + 150 * (2 * c[2]) -
                   324 * c[1]^2 * (2 * c[2]) -
                   1170 * c[1] * (3 * c[2]^2) -
                   1665 * (4 * c[2]^3))) -
    c[3] * ((8 - 6 * c[1]^4 - 2 * c[1]^6 + 144 * c[1] * c[2] -
               144 * c[1]^3 * c[2] - 108 * c[1]^5 * c[2] +
               720 * c[2]^2 - 540 * c[1]^2 * c[2]^2 -
               2178 * c[1]^4 * c[2]^2 + 2160 * c[1] * c[2]^3 -
               20952 * c[1]^3 * c[2]^3 + 9450 * c[2]^4 -
               106110 * c[1]^2 * c[2]^4 - 283500 * c[1] * c[2]^5 -
               330750 * c[2]^6)^((1/2) - 1) *
              ((1/2) * (144 * c[1] - 144 * c[1]^3 - 108 * c[1]^5 +
                          720 * (2 * c[2]) -
                          540 * c[1]^2 * (2 * c[2]) -
                          2178 * c[1]^4 * (2 * c[2]) +
                          2160 * c[1] * (3 * c[2]^2) -
                          20952 * c[1]^3 * (3 * c[2]^2) +
                          9450 * (4 * c[2]^3) -
                          106110 * c[1]^2 * (4 * c[2]^3) -
                          283500 * c[1] * (5 * c[2]^4) -
                          330750 * (6 * c[2]^5))))
  return(F)
}
