#' @title Calculate Theoretical Cumulative Probability
#'
#' @description This function calculates a cumulative probability using the theoretical power method cdf
#'     \eqn{F_p(Z)(p(z)) = F_p(Z)(p(z), F_Z(z))} up to \eqn{sigma * y + mu = delta}, where \eqn{y = p(z)}, after using
#'     \code{\link[SimMultiCorrData]{pdf_check}}.  If the given constants do not produce a valid power method pdf, a warning is given.
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param method the method used to find the constants.  "Fleishman" uses a third-order polynomial transformation and
#'     "Polynomial" uses Headrick's fifth-order transformation.
#' @param delta the value \eqn{sigma * y + mu}, where \eqn{y = p(z)}, at which to evaluate the cumulative probability
#' @param mu mean for the continuous variable
#' @param sigma standard deviation for the continuous variable
#' @param lower lower bound for integration of the standard normal variable Z that generates the continuous variable
#' @param upper upper bound for integration
#' @import stats
#' @import utils
#' @export
#' @keywords theoretical, statistics, cumulative, probability
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{pdf_check}}
#' @return A list with components:
#' @return \code{cumulative probability} the theoretical cumulative probability up to delta
#' @return \code{roots} the roots z that make \eqn{sigma * p(z) + mu = delta}
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
#'
#' @examples \dontrun{
#' # Beta(a = 4, b = 2) Distribution:
#' con <- find_constants(method = "Polynomial", skews = -0.467707, skurts = -0.375,
#'                       fifths = 1.403122, sixths = -0.426136)$constants
#' cdf_prob(c = con, method = "Polynomial", delta = 0.5)
#' }
cdf_prob <- function(c, method = c("Fleishman", "Polynomial"), delta = 0.5,
                     mu = 0, sigma = 1, lower = -1e06, upper = 1e06) {
  check <- suppressWarnings(pdf_check(c, method))
  if (check$valid.pdf == FALSE) {
    warning('This is NOT a valid pdf.')
  }
  cdf_root <- function(z, c, method, delta, mu, sigma) {
    if (method == "Fleishman") {
      y <- c[1] + c[2] * z + c[3] * z^2 + c[4] * z^3
    }
    if (method == "Polynomial") {
      y <- c[1] + c[2] * z + c[3] * z^2 + c[4] * z^3 + c[5] * z^4 +
        c[6] * z^5
    }
    return((sigma * y + mu) - delta)
  }
  R <- uniroot(cdf_root, lower = lower, upper = upper, c = c,
               method = method, delta = delta, mu = mu, sigma = sigma)
  cum_prob <- pnorm(R$root[1])
  return(list(cumulative_prob = cum_prob, roots = R$root))
}
