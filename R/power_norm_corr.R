#' @title Calculate Power Method Correlation
#'
#' @description This function calculates the correlation between a continuous variable, Y1, generated using a third or fifth-
#'     order polynomial transformation and the generating standard normal variable, Z1.  The power method correlation
#'     (described in Headrick & Kowalchuk, 2007) is given by:
#'     \eqn{\rho_{y1,z1} = c1 + 3*c3 + 15* c5}, where c5 = 0 if \code{method} = "Fleishman".  A value <= 0 indicates an invalid
#'     pdf and the signs of c1 and c3 should be reversed, which could still yield an invalid pdf.  All constants should
#'     be checked using \code{\link[SimMultiCorrData]{pdf_check}} to see if they generate a valid pdf.
#'
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#'
#' @param method the method used to find the constants.  "Fleishman" uses a third-order polynomial transformation and
#'     "Polynomial" uses Headrick's fifth-order transformation.
#'
#' @keywords correlation, Fleishman, Headrick
#'
#' @seealso \code{\link[SimMultiCorrData]{fleish}}, \code{\link[SimMultiCorrData]{poly}},
#'     \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{pdf_check}}
#'
#' @return A scalar equal to the correlation.
#' @references Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#' Non-normal Distributions. Computational Statistics & Data Analysis 40(4):685-711
#' (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35.
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3, 65-71.
#' @examples
#' # Beta(a = 4, b = 2) Distribution
#' power_norm_corr(c = c(0.108304, 1.104252, -0.123347, -0.045284, 0.005014,
#'                       0.001285),
#'                 method = "Polynomial")
#'
#' # Switch signs on c1 and c3 to get negative correlation (invalid pdf):
#' power_norm_corr(c = c(0.108304, -1.104252, -0.123347, 0.045284, 0.005014,
#'                       0.001285),
#'                 method = "Polynomial")
#'
#' @export
power_norm_corr <- function(c, method) {
  c <- as.numeric(c)
  if (method == "Fleishman") return(c[2] + 3 * c[4])
  if (method == "Polynomial") return(c[2] + 3 * c[4] + 15 * c[6])
}

