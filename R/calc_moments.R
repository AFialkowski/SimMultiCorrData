#' @title Find Standardized Cumulants of Data by Method of Moments
#'
#' @description This function uses the method of moments to calculate the mean, standard deviation, skewness,
#'     standardized kurtosis, and standardized fifth and sixth cumulants given a vector of data.  The result can be used
#'     as input to \code{\link[SimMultiCorrData]{find_constants}} or for data simulation.
#' @param x a vector of data
#' @export
#' @keywords cumulants, method of moments
#' @seealso \code{\link[SimMultiCorrData]{calc_fisherk}}, \code{\link[SimMultiCorrData]{calc_theory}},
#'          \code{\link[SimMultiCorrData]{find_constants}}
#' @return A vector of the mean, standard deviation, skewness, standardized kurtosis, and standardized fifth and sixth cumulants
#' @references
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis, 40(4):685-711. \doi{10.1016/S0167-9473(02)00072-5}.
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249. \doi{10.1080/10629360600605065}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' Kendall M & Stuart A (1977). The Advanced Theory of Statistics, 4th Edition. Macmillan, New York.
#'
#' @examples
#' x <- rgamma(n = 10000, 10, 10)
#' calc_moments(x)
calc_moments <- function(x) {
  # central moments
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m3 <- mean((x - m)^3)
  m4 <- mean((x - m)^4)
  m5 <- mean((x - m)^5)
  m6 <- mean((x - m)^6)
  s <- sqrt(m2)
  # central cumulants
  g1 <- m3/(s^3)
  g2 <- m4/(s^4) - 3
  g3 <- m5/(s^5) - 10 * g1
  g4 <- m6/(s^6) - 15 * g2 - 10 * g1^2 - 15
  stcums <- c(m, s, g1, g2, g3, g4)
  names(stcums) <- c("mean", "sd", "skew", "kurtosis", "fifth", "sixth")
  return(stcums)
}
