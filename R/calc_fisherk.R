#' @title Find Standardized Cumulants of Data based on Fisher's k-statistics
#'
#' @description This function uses Fisher's k-statistics to calculate the mean, standard deviation, skewness,
#'     standardized kurtosis, and standardized fifth and sixth cumulants given a vector of data.  The result can be used
#'     as input to \code{\link[SimMultiCorrData]{find_constants}} or for data simulation.
#' @param x a vector of data
#' @import stats
#' @import utils
#' @export
#' @keywords cumulants, Fisher
#' @seealso \code{\link[SimMultiCorrData]{calc_theory}}, \code{\link[SimMultiCorrData]{calc_moments}},
#'          \code{\link[SimMultiCorrData]{find_constants}}
#' @return A vector of the mean, standard deviation, skewness, standardized kurtosis, and standardized fifth and sixth cumulants
#' @references Fisher RA (1928). Moments and Product Moments of Sampling Distributions. Proc. London Math. Soc. 30, 199-238. \doi{10.1112/plms/s2-30.1.199}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}
#'
#' @examples
#' x <- rgamma(n = 10000, 10, 10)
#' calc_fisherk(x)
calc_fisherk <- function(x) {
  m <- mean(x)
  s <- sd(x)
  st_x <- scale(x, center = TRUE, scale = TRUE)
  n <- length(x)
  s2 <- n - 1
  s3 <- sum(st_x^3)
  s4 <- sum(st_x^4)
  s5 <- sum(st_x^5)
  s6 <- sum(st_x^6)
  g1 <- (n^2 * s3) / (n * (n - 1) * (n - 2))
  g2 <- ((n^3 + n^2) * s4 - 3 * (n^2 - n) *
            s2^2) / (n * (n - 1) * (n - 2) * (n - 3))
  g3 <- ((n^4 + 5 * n^3) * s5 - 10 * (n^3 - n^2) *
           s3 * s2) / (n * (n - 1) * (n - 2) * (n - 3) * (n - 4))
  g4 <- ((n^5 + 16 * n^4 + 11 * n^3 - 4 * n^2) * s6 - 15 *
           n * (n - 1)^2 * (n + 4) * s4 * s2 - 10 *
           (n^4 - 2 * n^3 + 5 * n^2 - 4 * n) * s3^2 + 30 *
           (n^3 - 3 * n^2 + 2 * n) * s2^3) / (n * (n - 1) * (n - 2) *
                                              (n - 3) * (n - 4) *
                                              (n - 5))
  stcums <- c(m, s, g1, g2, g3, g4)
  names(stcums) <- c("mean", "sd", "skew", "kurtosis", "fifth", "sixth")
  return(stcums)
}
