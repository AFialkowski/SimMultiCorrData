#' @title Calculate Theoretical Statistics for a Valid Power Method PDF
#'
#' @description This function calculates the 100*\code{alpha} percent symmetric trimmed mean (0 < \code{alpha} < 0.50), median,
#'     mode, and maximum height of a valid power method pdf, after using \code{\link[SimMultiCorrData]{pdf_check}}.  It will stop with
#'     an error if the pdf is invalid.
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param method the method used to find the constants.  "Fleishman" uses a third-order polynomial transformation and
#'     "Polynomial" uses Headrick's fifth-order transformation.
#' @param alpha proportion to be trimmed from the lower and upper ends of the power method pdf (default = 0.025)
#' @param mu mean for the continuous variable (default = 0)
#' @param sigma standard deviation for the continuous variable (default = 1)
#' @param lower lower bound for integration of the standard normal variable Z that generates the continuous variable (default = -10)
#' @param upper upper bound for integration (default = 10)
#' @param sub the number of subdivisions to use in the integration; if no result, try increasing sub (requires longer
#'     computation time; default = 1000)
#' @import stats
#' @import utils
#' @export
#' @keywords theoretical, statistics
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{pdf_check}}
#' @return A vector with components:
#' @return \code{trimmed_mean} the trimmed mean value
#' @return \code{median} the median value
#' @return \code{mode} the mode value
#' @return \code{max_height} the maximum pdf height
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
#' con <- find_constants(method = "Polynomial", skews = -0.467707,
#'                       skurts = -0.375, fifths = 1.403122,
#'                       sixths = -0.426136)$constants
#' stats_pdf(c = con, method = "Polynomial", alpha = 0.025)
#' }
stats_pdf <- function(c, method = c("Fleishman", "Polynomial"), alpha = 0.025,
                      mu = 0, sigma = 1, lower = -10, upper = 10, sub = 1000) {
  c <- as.numeric(c)
  check <- suppressWarnings(pdf_check(c, method))
  if (check$valid.pdf == FALSE) stop("This is NOT a valid pdf.")
  if (alpha <= 0 | alpha >= 0.5) stop("alpha value is not valid.")
  fy <- function(z, c, method) {
    if (method == "Fleishman") {
      return(-((1/sqrt(2 * pi)) * (exp(-0.5 * z^2) * (0.5 * (2 * z)))/(c[2] +
               2 * c[3] * z + 3 * c[4] * z^2) + ((1/sqrt(2 * pi)) * exp(-0.5 *
               z^2)) * (2 * c[3] + 3 * c[4] * (2 * z))/(c[2] + 2 * c[3] * z +
               3 * c[4] * z^2)^2))
    }
    if (method == "Polynomial") {
      return(-((1/sqrt(2 * pi)) * (exp(-0.5 * z^2) * (0.5 * (2 * z)))/(c[2] +
               2 * c[3] * z + 3 * c[4] * z^2 + 4 * c[5] * z^3 + 5 * c[6] *
               z^4) + ((1/sqrt(2 * pi)) * exp(-0.5 * z^2)) * (2 * c[3] +
               3 * c[4] * (2 * z) + 4 * c[5] * (3 * z^2) + 5 * c[6] *
              (4 * z^3))/(c[2] + 2 * c[3] * z + 3 * c[4] * z^2 + 4 * c[5] *
                            z^3 + 5 * c[6] * z^4)^2))
    }
  }
  f_max <- uniroot(fy, interval = c(lower, upper), c = c, method = method)
  z_tilde <- f_max$root
  if (method == "Fleishman") {
    mode <- sigma * (c[1] + c[2] * z_tilde + c[3] * z_tilde^2 +
                       c[4] * z_tilde^3) + mu
    max <- ((1/sqrt(2 * pi)) * exp(-0.5 * z_tilde^2))/(c[2] +
            2 * c[3] * z_tilde + 3 * c[4] * z_tilde^2)
  }
  if (method == "Polynomial") {
    mode <- sigma * (c[1] + c[2] * z_tilde + c[3] * z_tilde^2 +
                       c[4] * z_tilde^3 + c[5] * z_tilde^4 +
                       c[6] * z_tilde^5) + mu
    max <- ((1/sqrt(2 * pi)) * exp(-0.5 * z_tilde^2))/(c[2] +
            2 * c[3] * z_tilde + 3 * c[4] * z_tilde^2 +
              4 * c[5] * z_tilde^3 + 5 * c[6] * z_tilde^4)
  }
  median <- sigma * c[1] + mu
  phi2 <- function(z) (1/sqrt(2 * pi)) * exp(-0.5 * z^2)
  zlower <- qnorm(alpha)
  zupper <- qnorm(1 - alpha)
  p0 <- integrate(function(z, FUN = phi2) FUN(z), zlower, zupper,
                  subdivisions = sub, stop.on.error = FALSE)$value
  p1 <- integrate(function(z, FUN = phi2) z * FUN(z), zlower, zupper,
                  subdivisions = sub, stop.on.error = FALSE)$value
  p2 <- integrate(function(z, FUN = phi2) z^2 * FUN(z), zlower, zupper,
                  subdivisions = sub, stop.on.error = FALSE)$value
  p3 <- integrate(function(z, FUN = phi2) z^3 * FUN(z), zlower, zupper,
                  subdivisions = sub, stop.on.error = FALSE)$value
  p4 <- integrate(function(z, FUN = phi2) z^4 * FUN(z), zlower, zupper,
                  subdivisions = sub, stop.on.error = FALSE)$value
  p5 <- integrate(function(z, FUN = phi2) z^5 * FUN(z), zlower, zupper,
                  subdivisions = sub, stop.on.error = FALSE)$value
  if (method == "Fleishman") {
    ptotal <- c[1] * p0 + c[2] * p1 + c[3] * p2 + c[4] * p3
  }
  if (method == "Polynomial") {
    ptotal <- c[1] * p0 + c[2] * p1 + c[3] * p2 + c[4] * p3 + c[5] * p4 +
      c[6] * p5
  }
  tmean <- sigma * (1/(1 - 2 * alpha)) * ptotal + mu
  result <- c(tmean, median, mode, max)
  names(result) <- c("trimmed_mean", "median", "mode", "max_height")
  return(result)
}
