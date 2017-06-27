#' @title Find Theoretical Standardized Cumulants for Continuous Distributions
#'
#' @description This function calculates the theoretical mean, standard deviation, skewness,
#'     standardized kurtosis, and standardized fifth and sixth cumulants given either a Distribution name (plus up to 3
#'     parameters) or a pdf (with specified lower and upper support bounds).  The result can be used as input to
#'     \code{\link[SimMultiCorrData]{find_constants}} or for data simulation.
#'
#'     \bold{Note:} Due to the nature of the integration involved in calculating the standardized cumulants, the results are
#'     approximations.  Greater accuracy can be achieved by increasing the number of subdivisions (\code{sub}) used in the integration
#'     process.  However, the user may need to round the cumulants (i.e. using \code{round(x, 8)}) before using them in other functions
#'     (i.e. \code{find_constants}, \code{calc_lower_skurt}, \code{nonnormvar1}, \code{rcorrvar}, \code{rcorrvar2}) in order to achieve
#'     the desired results.  For example, in order to ensure that skew is exactly 0 for symmetric distributions.
#' @param Dist name of the distribution. The possible values are: "Beta", "Chisq", "Exponential", "F", "Gamma", "Gaussian",
#'     "Laplace", "Logistic", "Lognormal", "Pareto", "Rayleigh", "t", "Triangular", "Uniform", "Weibull".
#'     Please refer to the documentation for each package (i.e. \code{\link[stats]{dgamma}})
#'     for information on appropriate parameter inputs.  The pareto (see \code{\link[VGAM]{dpareto}}), generalized
#'     rayleigh (see \code{\link[VGAM]{dgenray}}), and laplace (see \code{\link[VGAM]{dlaplace}}) distributions
#'     come from the \code{\link[VGAM]{VGAM}} package.  The triangular (see \code{\link[triangle]{dtriangle}}) distribution
#'     comes from the \code{\link[triangle]{triangle}} package.
#' @param params a vector of parameters (up to 3) for the desired distribution (keep NULL if \code{fx} supplied instead)
#' @param fx a pdf input as a function of x only, i.e. fx <- function(x) 0.5*(x-1)^2; must return a scalar
#'     (keep NULL if Dist supplied instead)
#' @param lower the lower support bound for a supplied fx, else keep NULL
#' @param upper the upper support bound for a supplied fx, else keep NULL
#' @param sub the number of subdivisions to use in the integration; if no result, try increasing sub (requires longer
#'     computation time)
#' @export
#' @import stats
#' @import utils
#' @importFrom VGAM dpareto rpareto dgenray rgenray dlaplace rlaplace
#' @importFrom triangle dtriangle rtriangle
#' @keywords cumulants, theoretical
#' @seealso \code{\link[SimMultiCorrData]{calc_fisherk}}, \code{\link[SimMultiCorrData]{calc_moments}},
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
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}
#'
#' Thomas W. Yee (2017). VGAM: Vector Generalized Linear and Additive Models. R package version 1.0-3.
#'     \url{https://CRAN.R-project.org/package=VGAM}
#'
#' Rob Carnell (2016). triangle: Provides the Standard Distribution Functions for the Triangle Distribution. R package
#'     version 0.10. \url{https://CRAN.R-project.org/package=triangle}
#'
#' @examples
#' options(scipen = 999)
#'
#' # Pareto Distribution: params = c(alpha = scale, theta = shape)
#' calc_theory(Dist = "Pareto", params = c(1, 10))
#'
#' # Generalized Rayleigh Distribution: params = c(alpha = scale, mu/sqrt(pi/2) = shape)
#' calc_theory(Dist = "Rayleigh", params = c(0.5, 1))
#'
#' # Laplace Distribution: params = c(location, scale)
#' calc_theory(Dist = "Laplace", params = c(0, 1))
#'
#' # Triangle Distribution: params = c(a, b)
#' calc_theory(Dist = "Triangular", params = c(0, 1))
calc_theory <- function(Dist = c("Beta", "Chisq", "Exponential", "F", "Gamma",
                                 "Gaussian", "Laplace", "Logistic",
                                 "Lognormal", "Pareto", "Rayleigh", "t",
                                 "Triangular", "Uniform", "Weibull"),
                        params = NULL, fx = NULL, lower = NULL, upper = NULL,
                        sub = 1000) {
  if (!is.null(fx)) {
    m <- integrate(function(x, FUN = fx) x * FUN(x), lower, upper,
                   subdivisions = sub, stop.on.error = FALSE)$value
    m2 <- integrate(function(x, FUN = fx) (x - m)^2 * FUN(x), lower, upper,
                    subdivisions = sub, stop.on.error = FALSE)$value
    m3 <- integrate(function(x, FUN = fx) (x - m)^3 * FUN(x), lower, upper,
                    subdivisions = sub, stop.on.error = FALSE)$value
    m4 <- integrate(function(x, FUN = fx) (x - m)^4 * FUN(x), lower, upper,
                    subdivisions = sub, stop.on.error = FALSE)$value
    m5 <- integrate(function(x, FUN = fx) (x - m)^5 * FUN(x), lower, upper,
                    subdivisions = sub, stop.on.error = FALSE)$value
    m6 <- integrate(function(x, FUN = fx) (x - m)^6 * FUN(x), lower, upper,
                    subdivisions = sub, stop.on.error = FALSE)$value
  }
  if (is.null(fx)) {
    D <- data.frame(Dist = c("Beta", "Chisq", "Exponential", "F", "Gamma",
                            "Gaussian", "Laplace", "Logistic", "Lognormal",
                            "Pareto", "Rayleigh", "t", "Triangular",
                            "Uniform", "Weibull"),
                   pdf = c("dbeta", "dchisq", "dexp", "df", "dgamma",
                           "dnorm", "dlaplace", "dlogis", "dlnorm",
                           "dpareto", "dgenray", "dt", "dtriangle",
                           "dunif", "dweibull"),
                   fx = c("rbeta", "rchisq", "rexp", "rf", "rgamma",
                          "rnorm", "rlaplace", "rlogis", "rlnorm",
                          "rpareto", "rgenray", "rt", "rtriangle",
                          "runif", "rweibull"),
                   Lower = as.numeric(c(0, 0, 0, 0, 0, -Inf, -Inf, -Inf,
                                        0, params[1], 0, -Inf, params[1],
                                        params[1], 0)),
                   Upper = as.numeric(c(1, Inf, Inf, Inf, Inf, Inf, Inf,
                                        Inf, Inf, Inf, Inf, Inf, params[2],
                                        params[2], Inf)))
    i <- match(Dist, D$Dist)
    p <- as.character(D$pdf[i])
    if (length(params) == 1) {
      fx1 <- function(x) get(p)(x, params[1])
      m <- integrate(function(x, FUN = fx1) x * FUN(x), D$Lower[i], D$Upper[i],
                     subdivisions = sub, stop.on.error = FALSE)$value
      m2 <- integrate(function(x, FUN = fx1) (x - m)^2 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m3 <- integrate(function(x, FUN = fx1) (x - m)^3 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m4 <- integrate(function(x, FUN = fx1) (x - m)^4 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m5 <- integrate(function(x, FUN = fx1) (x - m)^5 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m6 <- integrate(function(x, FUN = fx1) (x - m)^6 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
    }
    if (length(params) == 2) {
      fx2 <- function(x) get(p)(x, params[1], params[2])
      m <- integrate(function(x, FUN = fx2) x * FUN(x), D$Lower[i], D$Upper[i],
                     subdivisions = sub, stop.on.error = FALSE)$value
      m2 <- integrate(function(x, FUN = fx2) (x - m)^2 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m3 <- integrate(function(x, FUN = fx2) (x - m)^3 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m4 <- integrate(function(x, FUN = fx2) (x - m)^4 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m5 <- integrate(function(x, FUN = fx2) (x - m)^5 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m6 <- integrate(function(x, FUN = fx2) (x - m)^6 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
    }
    if (length(params) == 3) {
      fx3 <- function(x) get(p)(x, params[1], params[2], params[3])
      m <- integrate(function(x, FUN = fx3) x * FUN(x), D$Lower[i], D$Upper[i],
                     subdivisions = sub, stop.on.error = FALSE)$value
      m2 <- integrate(function(x, FUN = fx3) (x - m)^2 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m3 <- integrate(function(x, FUN = fx3) (x - m)^3 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m4 <- integrate(function(x, FUN = fx3) (x - m)^4 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m5 <- integrate(function(x, FUN = fx3) (x - m)^5 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
      m6 <- integrate(function(x, FUN = fx3) (x - m)^6 * FUN(x), D$Lower[i],
                      D$Upper[i], subdivisions = sub,
                      stop.on.error = FALSE)$value
    }
  }
  s <- sqrt(m2)
  g1 <- m3/(s^3)
  g2 <- m4/(s^4) - 3
  g3 <- m5/(s^5) - 10 * g1
  g4 <- m6/(s^6) - 15 * g2 - 10 * g1^2 - 15
  stcums <- c(m, s, g1, g2, g3, g4)
  names(stcums) <- c("mean", "sd", "skew", "kurtosis", "fifth", "sixth")
  return(stcums)
}
