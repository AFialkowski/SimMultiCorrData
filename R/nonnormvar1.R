#' @title Generation of One Non-Normal Continuous Variable Using the Power Method
#'
#' @description This function simulates one non-normal continuous variable using either Fleishman's third-order or Headrick's
#'     fifth-order polynomial transformation.  If only one variable is desired and that variable is continuous, this function should
#'     be used.  Headrick & Kowalchuk (2007) outlined a general method for comparing a simulated distribution \eqn{Y} to a given
#'     theoretical distribution \eqn{Y^*}.  These steps can be found in the \bold{Comparison of Simulated Distribution to Theoretical
#'     Distribution or Empirical Data} vignette.
#'
#' @section Overview of Simulation Process:
#'     1) The constants are calculated for the continuous variable using \code{\link[SimMultiCorrData]{find_constants}}.  If no
#'     solutions are found that generate a valid power method pdf, the function will return constants that produce an invalid pdf
#'     (or a stop error if no solutions can be found).  Possible solutions include: 1) changing the seed, or 2) using a \code{Six} vector
#'     of sixth cumulant correction values (if \code{method} = "Polynomial").  Errors regarding constant
#'     calculation are the most probable cause of function failure.
#'
#'     2) An intermediate standard normal variate X of length n is generated.
#'
#'     3) Summary statistics are calculated.
#'
#' @param method the method used to generate the continuous variable.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param means mean for the continuous variable (default = 0)
#' @param vars variance (default = 1)
#' @param skews skewness value (default = 0)
#' @param skurts standardized kurtosis (kurtosis - 3, so that normal variables have a value of 0; default = 0)
#' @param fifths standardized fifth cumulant (not necessary for \code{method} = "Fleishman"; default = 0)
#' @param sixths standardized sixth cumulant (not necessary for \code{method} = "Fleishman"; default = 0)
#' @param Six a vector of correction values to add to the sixth cumulant if no valid pdf constants are found,
#'     ex: \code{Six = seq(0.01, 2, by = 0.01)}; if no correction is desired, set \code{Six = NULL} (default)
#' @param cstart initial values for root-solving algorithm (see \code{\link[BB]{multiStart}} for \code{method} = "Fleishman"
#'     or \code{\link[nleqslv]{nleqslv}} for \code{method} = "Polynomial").  If user specified,
#'     must be input as a matrix. If NULL and all 4 standardized cumulants (rounded to 3 digits) are within
#'     0.01 of those in Headrick's common distribution table (see \code{\link[SimMultiCorrData]{Headrick.dist}}
#'     data), uses his constants as starting values; else, generates n sets of random starting values from
#'     uniform distributions.
#' @param n the sample size (i.e. the length of the simulated variable; default = 10000)
#' @param seed the seed value for random number generation (default = 1234)
#' @importFrom psych describe
#' @import stats
#' @import utils
#' @import BB
#' @import nleqslv
#' @export
#' @keywords simulation, continuous, univariate, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}
#' @return A list with the following components:
#' @return \code{constants} a data.frame of the constants
#' @return \code{continuous_variable} a data.frame of the generated continuous variable
#' @return \code{summary_continuous} a data.frame containing a summary of the variable
#' @return \code{sixth_correction} the sixth cumulant correction value
#' @return \code{valid.pdf} "TRUE" if constants generate a valid pdf, else "FALSE"
#' @return \code{Constants_Time} the time in minutes required to calculate the constants
#' @return \code{Simulation_Time} the total simulation time in minutes
#' @references Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis 40(4):685-711
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3, 65-71.
#'
#' @examples \dontrun{
#' # Use Headrick & Kowalchuk's (2007) steps to compare a simulated exponential
#' # (mean = 2) variable to the theoretical exponential(mean = 2) density:
#'
#' # 1) Obtain the standardized cumulants:
#' stcums <- calc_theory(Dist = "Exponential", params = 0.5) # rate = 1/mean
#'
#' # 2) Simulate the variable:
#' H_exp <- nonnormvar1("Polynomial", means = 2, vars = 2, skews = stcums[3],
#'                     skurts = stcums[4], fifths = stcums[5],
#'                     sixths = stcums[6], Six = NULL, cstart = NULL,
#'                     n = 10000, seed = 1234)
#'
#' names(H_exp)
#' H_exp$constants
#' #           c0        c1       c2         c3          c4           c5
#' # 1 -0.3077396 0.8005605 0.318764 0.03350012 -0.00367481 0.0001587076
#'
#' # 3) Determine whether the constants produce a valid power method pdf:
#'
#' H_exp$valid.pdf
#' # [1] "TRUE"
#'
#' # 4) Select a critical value:
#'
#' # Let alpha = 0.05.
#' y_star <- qexp(1 - 0.05, rate = 0.5) # note that rate = 1/mean
#' y_star
#' # [1] 5.991465
#'
#' # 5) Solve m_{2}^{1/2} * p(z') + m_{1} - y* = 0 for z', where m_{1} and
#' # m_{2} are the 1st and 2nd moments of Y*:
#'
#' # The exponential(2) distribution has a mean and standard deviation equal
#' # to 2.
#' # Solve 2 * p(z') + 2 - y_star = 0 for z'
#' # p(z') = c0 + c1 * z' + c2 * z'^2 + c3 * z'^3 + c4 * z'^4 + c5 * z'^5
#'
#' f_exp <- function(z, c, y) {
#'   return(2 * (c[1] + c[2] * z + c[3] * z^2 + c[4] * z^3 + c[5] * z^4 +
#'               c[6] * z^5) + 2 - y)
#' }
#'
#' z_prime <- uniroot(f_exp, interval = c(-1e06, 1e06),
#'                    c = as.numeric(H_exp$constants), y = y_star)$root
#' z_prime
#' # [1] 1.644926
#'
#' # 6) Calculate 1 - Phi(z'), the corresponding probability for the
#' # approximation Y to Y* (i.e. 1 - Phi(z') = 0.05), and compare to target
#' # value alpha:
#'
#' 1 - pnorm(z_prime)
#' # [1] 0.04999249
#'
#' # 7) Plot a parametric graph of Y* and Y:
#'
#' plot_sim_pdf_theory(sim_y = as.numeric(H_exp$continuous_variable[, 1]),
#'                     overlay = TRUE, Dist = "Exponential", params = 0.5)
#'
#' # Note we can also plot the empirical cdf and show the cumulative
#' # probability up to y_star:
#'
#' plot_sim_cdf(sim_y = as.numeric(H_exp$continuous_variable[, 1]),
#'              calc_cprob = TRUE, delta = y_star)
#'
#' }

nonnormvar1 <- function(method = c("Fleishman", "Polynomial"), means = 0,
                        vars = 1, skews = 0, skurts = 0,
                        fifths = 0, sixths = 0, Six = NULL,
                        cstart = NULL, n = 10000, seed = 1234) {
  start.time <- Sys.time()
  start.time.constants <- Sys.time()
  cons <- suppressWarnings(find_constants(method = method, skews = skews,
                                          skurts = skurts, fifths = fifths,
                                          sixths = sixths, Six = Six,
                                          cstart = cstart, n = 25,
                                          seed = seed))
  if (length(cons) == 1)
    stop(paste("\nConstants can not be found for continuous variable 1.",
               sep = ""))
  con_solution <- cons$constants
  SixCorr <- cons$SixCorr1
  constants <- matrix(con_solution, nrow = 1)
  cat("Constants: Distribution ", 1, " \n")
  stop.time.constants <- Sys.time()
  set.seed(seed)
  X_cont <- matrix(scale(rnorm(n, 0, 1)), nrow = n, ncol = 1)
  Y <- matrix(1, nrow = n, ncol = 1)
  Yb <- matrix(1, nrow = n, ncol = 1)
  if (method == "Fleishman") {
    Y[, 1] <- constants[1, 1] + constants[1, 2] * X_cont[, 1] +
      constants[1, 3] * X_cont[, 1]^2 + constants[1, 4] * X_cont[, 1]^3
    colnames(constants) <- c("c0", "c1", "c2", "c3")
  }
  if (method == "Polynomial") {
    Y[, 1] <- constants[1, 1] + constants[1, 2] * X_cont[, 1] +
      constants[1, 3] * X_cont[, 1]^2 + constants[1, 4] * X_cont[, 1]^3 +
      constants[1, 5] * X_cont[, 1]^4 + constants[1, 6] * X_cont[, 1]^5
    colnames(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5")
  }
  Yb[, 1] <- means + sqrt(vars) * Y[, 1]
  cont_sum <- describe(Yb, type = 1)
  sim_fifths <- matrix(calc_moments(Yb[, 1])[5], nrow = 1, ncol = 1)
  sim_sixths <- matrix(calc_moments(Yb[, 1])[6], nrow = 1, ncol = 1)
  cont_sum <- as.data.frame(cbind(1, cont_sum[, -c(1, 6, 7, 10, 13)],
                                  sim_fifths, sim_sixths))
  colnames(cont_sum) <- c("Distribution", "n", "mean", "sd", "median", "min",
                          "max", "skew", "skurtosis", "fifth", "sixth")
  stop.time <- Sys.time()
  Time.constants <- round(difftime(stop.time.constants, start.time.constants,
                                   units = "min"), 3)
  cat("\nConstants calculation time:", Time.constants, "minutes \n")
  Time <- round(difftime(stop.time, start.time, units = "min"), 3)
  cat("Total Simulation time:", Time, "minutes \n")
  result <- list(constants = as.data.frame(constants),
                 continuous_variable = as.data.frame(Yb),
                 summary_continuous = cont_sum, sixth_correction = SixCorr,
                 valid.pdf = cons$valid, Constants_Time = Time.constants,
                 Simulation_Time = Time)
  result
}
