#' @title Find Lower Boundary of Standardized Kurtosis for Polynomial Transformation
#'
#' @description This function calculates the lower boundary of standardized kurtosis for Fleishman's Third-Order (\code{method} = "Fleishman",
#'     \doi{10.1007/BF02293811}) or Headrick's Fifth-Order (\code{method} = "Polynomial", \doi{10.1016/S0167-9473(02)00072-5}), given values of skewness and standardized fifth and sixth
#'     cumulants.  It uses \code{\link[nleqslv]{nleqslv}} to search for solutions to the multi-constraint Lagrangean expression in either
#'     \code{\link[SimMultiCorrData]{fleish_skurt_check}} or \code{\link[SimMultiCorrData]{poly_skurt_check}}.  When Headrick's method
#'     is used (\code{method} = "Polynomial"), if no solutions converge and a vector of sixth cumulant correction values (\code{Six}) is
#'     provided, the smallest value is found that yields solutions.  Otherwise, the function stops with an error.
#'
#'     Each set of constants is checked for a positive correlation with the underlying normal variable
#'     (using \code{\link[SimMultiCorrData]{power_norm_corr}}) and a valid power method pdf (using \code{\link[SimMultiCorrData]{pdf_check}}).
#'     If the correlation is <= 0, the signs of c1 and c3 are reversed (for \code{method} = "Fleishman"),
#'     or c1, c3, and c5 (for \code{method} = "Polynomial").  It will return a kurtosis value with constants that yield in invalid pdf
#'     if no other solutions can be found (\code{valid.pdf} = "FALSE").  If a vector of kurtosis correction values (\code{Skurt}) is provided, the function
#'     finds the smallest value that produces a kurtosis with constants that yield a valid pdf.  If valid pdf constants
#'     still can not be found, the original invalid pdf constants (calculated without a correction) will be provided.  If no solutions
#'     can be found, an error is given and the function stops.  Please note that this function can take
#'     considerable computation time, depending on the number of starting values (n) and lengths of kurtosis (\code{Skurt}) and sixth cumulant
#'     (\code{Six}) correction vectors.  Different seeds should be tested to see if a lower boundary can be found.
#'
#' @section Notes on Fleishman Method:
#'     The Fleishman method \emph{can not generate valid power method distributions} with a ratio of \eqn{skew^2/skurtosis > 9/14}, where skurtosis is kurtosis - 3.
#'     This prevents the method from being used for any of the \emph{Chi-squared distributions}, which have a constant ratio of
#'     \eqn{skew^2/skurtosis = 2/3}.
#'
#'     \bold{Symmetric Distributions:} All symmetric distributions (which have skew = 0) possess the same lower kurtosis boundary.  This is solved for
#'     using \code{\link[stats]{optimize}} and the equations in Headrick & Sawilowsky (2002, \doi{10.3102/10769986025004417}).
#'     The result will always be: c0 = 0, c1 = 1.341159,
#'     c2 = 0, c3 = -0.1314796, and minimum standardized kurtosis = -1.151323.  Note that this set of constants does NOT generate a valid
#'     power method pdf.  If a \code{Skurt} vector of kurtosis correction values is provided, the function will find the smallest addition that yields a
#'     valid pdf.  This value is 1.16, giving a lower kurtosis boundary of 0.008676821.
#'
#'     \bold{Asymmetric Distributions:} Due to the square roots involved in the calculation of the lower kurtosis boundary (see Headrick & Sawilowsky, 2002),
#'     this function uses the absolute value of the skewness.  If the true skewness is less than zero, the signs on the constants c0 and c2 are
#'     switched after calculations (which changes skewness from positive to negative without affecting kurtosis).
#'
#'     \bold{Verification of Minimum Kurtosis:} Since differentiability is a local property, it is possible to obtain a local, instead of a global, minimum.
#'     For the Fleishman method, Headrick & Sawilowsky (2002) explain that since the equation for kurtosis is not "quasiconvex on the domain
#'     consisting only of the nonnegative orthant," second-order conditions must be verified.  The solutions for
#'     lambda, c1, and c3 generate a global kurtosis minimum if and only if the determinant of a bordered Hessian is less than zero.  Therefore,
#'     this function first obtains the solutions to the Lagrangean expression in \code{\link[SimMultiCorrData]{fleish_skurt_check}} for a given
#'     skewness value.  These are used to calculate the standardized kurtosis, the constants c1 and c3, and the Hessian determinant
#'     (using \code{\link[SimMultiCorrData]{fleish_Hessian}}).  If this determinant is less than zero, the kurtosis is indicated as a minimum.
#'     The constants c0, c1, c2, and c3 are checked to see if they yield a continuous variable with a positive correlation with the generating
#'     standard normal variable (using \code{\link[SimMultiCorrData]{power_norm_corr}}).  If not, the signs of c1 and c3 are switched.
#'     The final set of constants is checked to see if they generate a valid power method pdf (using \code{\link[SimMultiCorrData]{pdf_check}}).
#'     If a \code{Skurt} vector of kurtosis correction values is provided, the function will find the smallest value that yields a
#'     valid pdf.
#'
#' @section Notes on Headrick's Method:
#'     The \emph{sixth cumulant correction vector} (\code{Six}) may be used in order to aid in obtaining solutions which converge.  The calculation methods
#'     are the same for symmetric or asymmetric distributions, and for positive or negative skew.
#'
#'     \bold{Verification of Minimum Kurtosis:} For the fifth-order approximation, Headrick (2002, \doi{10.1016/S0167-9473(02)00072-5}) states
#'     "it is assumed that the hypersurface of the objective function [for the kurtosis
#'     equation] has the appropriate (quasiconvex) configuration."  This assumption alleviates the need to check second-order conditions.
#'     Headrick discusses steps he took to verify the kurtosis solution was in fact a minimum, including: 1) substituting the constant solutions
#'     back into the 1st four Lagrangean constraints to ensure the results are zero, 2) substituting the skewness, kurtosis solution, and
#'     standardized fifth and sixth cumulants back into the fifth-order equations to ensure the same constants are produced
#'     (i.e. using \code{\link[SimMultiCorrData]{find_constants}}), and 3) searching for values below the kurtosis solution that solve the
#'     Lagrangean equation.  This function ensures steps 1 and 2 by the nature of the root-solving algorithm of \code{\link[nleqslv]{nleqslv}}.
#'     Using a sufficiently large n (and, if necessary, executing the function for different seeds) makes step 3 unnecessary.
#'
#' @section Reasons for Function Errors:
#'     The most likely cause for function errors is that no solutions to \code{\link[SimMultiCorrData]{fleish_skurt_check}} or
#'     \code{\link[SimMultiCorrData]{poly_skurt_check}} converged.  If this happens,
#'     the simulation will stop.  Possible solutions include: a) increasing the number of initial starting values (\code{n}),
#'     b) using a different seed, or c) specifying a \code{Six} vector of sixth cumulant correction values (for \code{method} = "Polynomial").
#'     If the standardized cumulants are obtained from \code{calc_theory}, the user may need to use rounded values as inputs (i.e.
#'     \code{skews = round(skews, 8)}).  Due to the nature of the integration involved in \code{calc_theory}, the results are
#'     approximations.  Greater accuracy can be achieved by increasing the number of subdivisions (\code{sub}) used in the integration
#'     process.  For example, in order to ensure that skew is exactly 0 for symmetric distributions.
#'
#' @param method the method used to find the constants.  "Fleishman" uses a third-order polynomial transformation and
#'     requires only a skewness input.  "Polynomial" uses Headrick's fifth-order transformation and requires skewness plus standardized
#'     fifth and sixth cumulants.
#' @param skews the skewness value
#' @param fifths the standardized fifth cumulant (if \code{method} = "Fleishman", keep NULL)
#' @param sixths the standardized sixth cumulant (if \code{method} = "Fleishman", keep NULL)
#' @param Skurt a vector of correction values to add to the lower kurtosis boundary if the constants yield an invalid pdf,
#'     ex: \code{Skurt} = seq(0.1, 10, by = 0.1)
#' @param Six a vector of correction values to add to the sixth cumulant if no solutions converged,
#'     ex: \code{Six} = seq(0.05, 2, by = 0.05)
#' @param xstart initial value for root-solving algorithm (see \code{\link[nleqslv]{nleqslv}}).  If user specified,
#'     must be input as a matrix. If NULL generates n sets of random starting values from
#'     uniform distributions.
#' @param seed the seed value for random starting value generation (default = 104)
#' @param n the number of initial starting values to use (default = 50).  More starting values require more calculation time.
#' @export
#' @import nleqslv
#' @import BB
#' @import stats
#' @import utils
#' @keywords kurtosis, boundary, Fleishman, Headrick
#' @seealso \code{\link[nleqslv]{nleqslv}}, \code{\link[SimMultiCorrData]{fleish_skurt_check}},
#'     \code{\link[SimMultiCorrData]{fleish_Hessian}}, \code{\link[SimMultiCorrData]{poly_skurt_check}},
#'     \code{\link[SimMultiCorrData]{power_norm_corr}}, \code{\link[SimMultiCorrData]{pdf_check}},
#'     \code{\link[SimMultiCorrData]{find_constants}}
#' @return A list with components:
#' @return \code{Min}    a data.frame containing the skewness, fifth and sixth standardized cumulants (if \code{method} = "Polynomial"), constants,
#'     a valid.pdf column indicating whether or not the constants generate a valid power method pdf, and the minimum value of standardized
#'     kurtosis ("skurtosis")
#' @return \code{C}    a data.frame of valid power method pdf solutions, containing the skewness, fifth and sixth standardized cumulants
#'     (if \code{method} = "Polynomial"), constants, a valid.pdf column indicating TRUE, and all values of standardized kurtosis ("skurtosis").
#'     If the Lagrangean equations yielded valid pdf solutions, this will also include the lambda values, and for \code{method} = "Fleishman", the
#'     Hessian determinant and a minimum column indicating TRUE if the solutions give a minimum kurtosis.  If the Lagrangean equations
#'     yielded invalid pdf solutions, this data.frame contains constants calculated from \code{\link[SimMultiCorrData]{find_constants}}
#'     using the kurtosis correction.
#' @return \code{Invalid.C}    if the Lagrangean equations yielded invalid pdf solutions, a data.frame containing the skewness, fifth and sixth
#'     standardized cumulants (if \code{method} = "Polynomial"), constants, lambda values, a valid.pdf column indicating FALSE, and all values
#'     of standardized kurtosis ("skurtosis").  If \code{method} = "Fleishman", also the
#'     Hessian determinant and a minimum column indicating TRUE if the solutions give a minimum kurtosis.
#' @return \code{Time}    the total calculation time in minutes
#' @return \code{start}    a matrix of starting values used in root-solver
#' @return \code{SixCorr1}    if Six is specified, the sixth cumulant correction required to achieve converged solutions
#' @return \code{SkurtCorr1}    if Skurt is specified, the kurtosis correction required to achieve a valid power method pdf (or the maximum value
#'     attempted if no valid pdf solutions could be found)
#' @references
#' Berend H (2017). nleqslv: Solve Systems of Nonlinear Equations. R package version 3.2.
#'     \url{https://CRAN.R-project.org/package=nleqslv}
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532. \doi{10.1007/BF02293811}.
#'
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis, 40(4):685-711. \doi{10.1016/S0167-9473(02)00072-5}.
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3(1), 65-71. \doi{10.22237/jmasm/1083370080}.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249. \doi{10.1080/10629360600605065}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35. \doi{10.1007/BF02294317}.
#'
#' Headrick TC, Sawilowsky SS (2002). Weighted Simplex Procedures for Determining Boundary Points and Constants for the
#'     Univariate and Multivariate Power Methods. Journal of Educational and Behavioral Statistics, 25, 417-436. \doi{10.3102/10769986025004417}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' @examples \dontrun{
#'
#' # This example takes considerable computation time.
#'
#' # Reproduce Headrick's Table 2 (2002, p.698): note the seed here is 104.
#' # If you use seed = 1234, you get higher Headrick kurtosis values for V7 and V9.
#' # This shows the importance of trying different seeds.
#'
#' options(scipen = 999)
#'
#' V1 <- c(0, 0, 28.5)
#' V2 <- c(0.24, -1, 11)
#' V3 <- c(0.48, -2, 6.25)
#' V4 <- c(0.72, -2.5, 2.5)
#' V5 <- c(0.96, -2.25, -0.25)
#' V6 <- c(1.20, -1.20, -3.08)
#' V7 <- c(1.44, 0.40, 6)
#' V8 <- c(1.68, 2.38, 6)
#' V9 <- c(1.92, 11, 195)
#' V10 <- c(2.16, 10, 37)
#' V11 <- c(2.40, 15, 200)
#'
#' G <- as.data.frame(rbind(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11))
#' colnames(G) <- c("g1", "g3", "g4")
#'
#' # kurtosis correction vector (used in case of invalid power method pdf constants)
#' Skurt <- seq(0.01, 2, 0.01)
#'
#' # sixth cumulant correction vector (used in case of no converged solutions for
#' # method = "Polynomial")
#' Six <- seq(0.1, 10, 0.1)
#'
#' # Fleishman's Third-order transformation
#' F_lower <- list()
#' for (i in 1:nrow(G)) {
#'   F_lower[[i]] <- calc_lower_skurt("Fleishman", G[i, 1], Skurt = Skurt,
#'                                    seed = 104)
#' }
#'
#' # Headrick's Fifth-order transformation
#' H_lower <- list()
#' for (i in 1:nrow(G)) {
#'   H_lower[[i]] <- calc_lower_skurt("Polynomial", G[i, 1], G[i, 2], G[i, 3],
#'                                    Skurt = Skurt, Six = Six, seed = 104)
#' }
#'
#' # Approximate boundary from PoisBinOrdNonNor
#' PBON_lower <- G$g1^2 - 2
#'
#' # Compare results:
#' # Note: 1) the lower Headrick kurtosis boundary for V4 is slightly lower than the
#' #          value found by Headrick (-0.480129), and
#' #       2) the approximate lower kurtosis boundaries used in PoisBinOrdNonNor are
#' #          much lower than the actual Fleishman boundaries, indicating that the
#' #          guideline is not accurate.
#' Lower <- matrix(1, nrow = nrow(G), ncol = 12)
#' colnames(Lower) <- c("skew", "fifth", "sixth", "H_valid.skurt",
#'                      "F_valid.skurt", "H_invalid.skurt", "F_invalid.skurt",
#'                      "PBON_skurt", "H_skurt_corr", "F_skurt_corr",
#'                      "H_time", "F_time")
#'
#' for (i in 1:nrow(G)) {
#'   Lower[i, 1:3] <- as.numeric(G[i, 1:3])
#'   Lower[i, 4] <- ifelse(H_lower[[i]]$Min[1, "valid.pdf"] == "TRUE",
#'                         H_lower[[i]]$Min[1, "skurtosis"], NA)
#'   Lower[i, 5] <- ifelse(F_lower[[i]]$Min[1, "valid.pdf"] == "TRUE",
#'                         F_lower[[i]]$Min[1, "skurtosis"], NA)
#'   Lower[i, 6] <- min(H_lower[[i]]$Invalid.C[, "skurtosis"])
#'   Lower[i, 7] <- min(F_lower[[i]]$Invalid.C[, "skurtosis"])
#'   Lower[i, 8:12] <- c(PBON_lower[i], H_lower[[i]]$SkurtCorr1,
#'                       F_lower[[i]]$SkurtCorr1,
#'                       H_lower[[i]]$Time, F_lower[[i]]$Time)
#' }
#' Lower <- as.data.frame(Lower)
#' print(Lower[, 1:8], digits = 4)
#'
#' #    skew fifth  sixth H_valid.skurt F_valid.skurt H_invalid.skurt F_invalid.skurt PBON_skurt
#' # 1  0.00  0.00  28.50       -1.0551      0.008677         -1.3851         -1.1513    -2.0000
#' # 2  0.24 -1.00  11.00       -0.8600      0.096715         -1.2100         -1.0533    -1.9424
#' # 3  0.48 -2.00   6.25       -0.5766      0.367177         -0.9266         -0.7728    -1.7696
#' # 4  0.72 -2.50   2.50       -0.1319      0.808779         -0.4819         -0.3212    -1.4816
#' # 5  0.96 -2.25  -0.25        0.4934      1.443567          0.1334          0.3036    -1.0784
#' # 6  1.20 -1.20  -3.08        1.2575      2.256908          0.9075          1.1069    -0.5600
#' # 7  1.44  0.40   6.00            NA      3.264374          1.7758          2.0944     0.0736
#' # 8  1.68  2.38   6.00            NA      4.452011          2.7624          3.2720     0.8224
#' # 9  1.92 11.00 195.00        5.7229      5.837442          4.1729          4.6474     1.6864
#' # 10 2.16 10.00  37.00            NA      7.411697          5.1993          6.2317     2.6656
#' # 11 2.40 15.00 200.00            NA      9.182819          6.6066          8.0428     3.7600
#'
#' Lower[, 9:12]
#'
#' #    H_skurt_corr F_skurt_corr H_time F_time
#' # 1          0.33         1.16  1.757  8.227
#' # 2          0.35         1.15  1.566  8.164
#' # 3          0.35         1.14  1.630  6.321
#' # 4          0.35         1.13  1.537  5.568
#' # 5          0.36         1.14  1.558  5.540
#' # 6          0.35         1.15  1.602  6.619
#' # 7          2.00         1.17  9.088  8.835
#' # 8          2.00         1.18  9.425 11.103
#' # 9          1.55         1.19  6.776 14.364
#' # 10         2.00         1.18 11.174 15.382
#' # 11         2.00         1.14 10.567 18.184
#'
#' # The 1st 3 columns give the skewness and standardized fifth and sixth cumulants.
#' # "H_valid.skurt" gives the lower kurtosis boundary that produces a valid power method pdf
#' #     using Headrick's approximation, with the kurtosis addition given in the "H_skurt_corr"
#' #     column if necessary.
#' # "F_valid.skurt" gives the lower kurtosis boundary that produces a valid power method pdf
#' #     using Fleishman's approximation, with the kurtosis addition given in the "F_skurt_corr"
#' #     column if necessary.
#' # "H_invalid.skurt" gives the lower kurtosis boundary that produces an invalid power method
#' #     pdf using Headrick's approximation, without the use of a kurtosis correction.
#' # "F_valid.skurt" gives the lower kurtosis boundary that produces an invalid power method
#' #     pdf using Fleishman's approximation, without the use of a kurtosis correction.
#' # "PBON_skurt" gives the lower kurtosis boundary approximation used in the PoisBinOrdNonNor
#' #     package.
#' # "H_time" gives the computation time (minutes) for Headrick's method.
#' # "F_time" gives the computation time (minutes) for Fleishman's method.
#'
#' }
calc_lower_skurt <- function(method = c("Fleishman", "Polynomial"),
                             skews = NULL, fifths = NULL, sixths = NULL,
                             Skurt = NULL, Six = NULL, xstart = NULL,
                             seed = 104, n = 50) {
  start.time <- Sys.time()
  if (method == "Fleishman") {
    SkurtCorr1 <- NULL
    error1 <- "Lower boundary could not be found.
    Try more starting values (increase n) or a different seed.
    Also check the cumulant values.\n"
    error2 <- "Only invalid power method constants could be found.
    Try a larger n, a different seed, or a kurtosis correction vector.
    Also check the cumulant values.\n"
    skews0 <- skews
    skews <- abs(skews)
    if (skews[1] == 0) {
      gamma2 <- function(c3) 24 * ((sqrt(1 - 6 * c3^2) - 3 * c3) *
                                     c3 + c3^2 * (12 + 48 *
                                               (sqrt(1 - 6 * c3^2) -
                                                  3 * c3) * c3 +
                                               225 * c3^2))
      c3 <- suppressWarnings(optimize(gamma2,
                                      interval = c(-1e06,
                                              sqrt(1/6)))$minimum)
      g2 <- gamma2(c3)
      c1 <- sqrt(1 - 6 * c3^2) - 3 * c3
      C <- matrix(1, nrow = 1, ncol = 6)
      C[1, ] <- c(0, 0, c1, 0, c3, g2)
      C <- data.frame(C, suppressWarnings(pdf_check(c(0, c1, 0, c3),
                                   method)$valid.pdf))
      colnames(C) <- c("skew", "c0", "c1", "c2", "c3", "skurtosis",
                       "valid.pdf")
      constants2 <- C[C$valid.pdf == "TRUE", , drop = FALSE]
      if (length(Skurt) == 0) {
        stop.time <- Sys.time()
        Time <- round(difftime(stop.time, start.time,
                               units = "min"), 3)
        return(list(Min = C, C = C, Invalid.C = C, Time = Time,
                    SkurtCorr1 = SkurtCorr1))
      } else {
        v <- 1
        while (nrow(constants2) == 0) {
          con_solution <-
            suppressWarnings(find_constants(method = method,
                                            skews = skews,
                                            skurts =
                                              min(C$skurtosis) +
                                              Skurt[v],
                                            fifths = NULL,
                                            sixths = NULL,
                                            Six = NULL,
                                            cstart = NULL,
                                            n = n,
                                            seed = seed))
          if (length(con_solution) != 1) {
            if (con_solution$valid == "TRUE") {
              constants2 <- matrix(con_solution$constants,
                                   nrow = 1, ncol =
                                     length(con_solution$constants))
            }
          }
          v <- v + 1
          if (v > length(Skurt)) break
        }
        SkurtCorr1 <- Skurt[v - 1]
        if (nrow(constants2) != 0) {
          constants2 <-  data.frame(skews, constants2, "TRUE",
                                    min(C$skurtosis) +
                                      Skurt[v - 1])
          colnames(constants2) <- c("skew", "c0", "c1",
                                    "c2", "c3", "valid.pdf",
                                    "skurtosis")
          stop.time <- Sys.time()
          Time <- round(difftime(stop.time, start.time,
                                 units = "min"), 3)
          return(list(Min = constants2, C = constants2,
                      Invalid.C = C, Time = Time,
                      SkurtCorr1 = SkurtCorr1))
        } else {
          stop.time <- Sys.time()
          Time <- round(difftime(stop.time, start.time,
                                 units = "min"), 3)
          return(list(Min = C, C = NULL, Invalid.C = C,
                      Time = Time, SkurtCorr1 = SkurtCorr1))
        }
      }
    }
    if (skews[1] != 0) {
      if (length(xstart) == 0) {
        set.seed(seed)
        cstart1 <- runif(n, min = 0.5, max = 1.35)
        cstart2 <- runif(n, min = -0.15, max = 0.05)
        lstart1 <- runif(n, min = 0, max = 10)
        xstart <- cbind(cstart1, cstart2, lstart1)
      }
      test <- numeric(nrow(xstart))
      for (i in 1:nrow(xstart)) {
        test[i] <- ifelse(!is.na(fleish_skurt_check(xstart[i, ],
                                                    a = skews)[1]),
                          TRUE, FALSE)
      }
      xstart2 <- xstart[which(test == TRUE), , drop = FALSE]
      if (is.null(xstart2)) {
        stop(error1)
      } else {
        converged <- NULL
        for (i in 1:nrow(xstart2)) {
          nl_solution <- nleqslv(x = xstart2[i, ],
                                 fn = fleish_skurt_check,
                                 a = skews, method = "Broyden",
                                 control = list(ftol = 1e-05))
          if (nl_solution$termcd == 1) {
            converged <- rbind(converged, nl_solution$x)
          }
        }
        if (is.null(converged)) {
          stop(error1)
        } else {
          converged <- converged[!duplicated(converged), , drop = FALSE]
          constants <- data.frame(converged,
                                  rep("FALSE", nrow(converged)),
                                  rep(0, nrow(converged)),
                                  rep(0, nrow(converged)),
                                  rep(0, nrow(converged)),
                                  rep(0, nrow(converged)))
          colnames(constants) <- c("c1", "c3", "lambda",
                                   "valid.pdf", "H_det",
                                   "skurtosis", "c2", "c0")
          constants$valid.pdf <- as.character(constants$valid.pdf)
          for (i in 1:nrow(constants)) {
            constants$H_det[i] <-
              fleish_Hessian(as.numeric(constants[i, 1:3]))$H_det
            constants$skurtosis[i] <-
              12 * (1 - constants[i, 1]^4 + 24 * constants[i, 1] *
                      constants[i, 2] - 34 * constants[i, 1]^3 *
                      constants[i, 2] + 150 * constants[i, 2]^2 -
                      324 * constants[i, 1]^2 * constants[i, 2]^2 -
                      1170 * constants[i, 1] * constants[i, 2]^3 -
                      1665 * constants[i, 2]^4)
            constants$c2[i] <-
              skews[1]/(2 * (constants[i, 1]^2 +
                               24 * constants[i, 1] *
                               constants[i, 2] +
                               105 * constants[i, 2]^2 + 2))
            constants$c0[i] <- -constants$c2[i]
            if (skews0 < 0) {
              constants$c2[i] <- -constants$c2[i]
              constants$c0[i] <- -constants$c0[i]
            }
          }
          constants <- data.frame(skew = rep(skews0,
                                             nrow(constants)),
                                  c0 = constants$c0,
                                  c1 = constants$c1,
                                  c2 = constants$c2,
                                  c3 = constants$c3,
                                  lambda = constants$lambda,
                                  valid.pdf = constants$valid.pdf,
                                  skurtosis = constants$skurtosis,
                                  H_det = constants$H_det)
          constants$Min <- ifelse(constants$H_det < 0,
                                  "yes", "no")
          for (i in 1:nrow(constants)) {
            if (power_norm_corr(as.numeric(constants[i, 2:5]),
                                method = method) < 0) {
              constants[i, 3] <- -constants[i, 3]
              constants[i, 5] <- -constants[i, 5]
            }
            check <-
              suppressWarnings(pdf_check(as.numeric(constants[i, 2:5]),
                                       method = method)$valid.pdf)
            if (check[1] == TRUE) {
              constants$valid.pdf[i] <- "TRUE"
            } else {
              constants$valid.pdf[i] <- "FALSE"
            }
          }
          constants2 <- constants[constants$valid.pdf == "TRUE", , drop = FALSE]
          if (nrow(constants2) != 0) {
            m <- min(constants2$skurtosis)
            C2 <- constants2[constants2$skurtosis == m, , drop = FALSE]
            stop.time <- Sys.time()
            Time <- round(difftime(stop.time, start.time,
                                   units = "min"), 3)
            return(list(Min = C2, C = constants2, Invalid.C = NULL,
                        start = xstart2, Time = Time,
                        SkurtCorr1 = SkurtCorr1))
          } else {
            if (length(Skurt) == 0) {
              m <- min(constants$skurtosis)
              C2 <- constants[constants$skurtosis == m, , drop = FALSE]
              stop.time <- Sys.time()
              Time <- round(difftime(stop.time, start.time,
                                     units = "min"), 3)
              cat(error2)
              return(list(Min = C2, C = NULL,
                          Invalid.C = constants, start = xstart2,
                          Time = Time, SkurtCorr1 = SkurtCorr1))
            } else {
              v <- 1
              while (nrow(constants2) == 0) {
                con_solution <-
                  suppressWarnings(find_constants(method = method,
                                                  skews = skews,
                                                  skurts =
                                                    min(constants$skurtosis) +
                                                    Skurt[v],
                                                  fifths = NULL,
                                                  sixths = NULL,
                                                  Six = NULL, cstart = NULL,
                                                  n = n, seed = seed))
                if (length(con_solution) != 1) {
                  if (con_solution$valid == "TRUE") {
                    constants2 <- matrix(con_solution$constants, nrow = 1,
                                         ncol = length(con_solution$constants))
                  }
                }
                v <- v + 1
                if (v > length(Skurt)) break
              }
              SkurtCorr1 <- Skurt[v - 1]
            }
            if (nrow(constants2) != 0) {
              constants2 <-  data.frame(skews, constants2, "TRUE",
                                      min(constants$skurtosis) + Skurt[v - 1])
              colnames(constants2) <- c("skew", "c0", "c1", "c2", "c3",
                                        "valid.pdf", "skurtosis")
              m <- min(constants2$skurtosis)
              C2 <- constants2[constants2$skurtosis == m, , drop = FALSE]
              stop.time <- Sys.time()
              Time <- round(difftime(stop.time, start.time,
                                     units = "min"), 3)
              return(list(Min = C2, C = constants2, Invalid.C = constants,
                          start = xstart2,
                          Time = Time, SkurtCorr1 = SkurtCorr1))
            }
            if (nrow(constants2) == 0) {
              m <- min(constants$skurtosis)
              C2 <- constants[constants$skurtosis == m, , drop = FALSE]
              stop.time <- Sys.time()
              Time <- round(difftime(stop.time, start.time,
                                     units = "min"), 3)
              cat(error2)
              return(list(Min = C2, C = NULL, Invalid.C = constants,
                          start = xstart2,
                          Time = Time, SkurtCorr1 = SkurtCorr1))
            }
          }
        }
      }
    }
  }
  if (method == "Polynomial") {
    error1 <- "Lower boundary could not be found.
               Try more starting values (increase n)
               or a different seed or Six vector.
               Also verify standardized cumulant values.\n"
    error2 <- "Only invalid pdf constants could be found.
               Try more starting values (increase n)
               or a different seed or Skurt vector.
               Also verify standardized cumulant values.\n"
    SixCorr1 <- NULL
    SkurtCorr1 <- NULL
    poly_skurt <- function(c) {
      c <- as.numeric(c)
      pskurtosis <- 24 * (2 * c[2]^4 + 96 * c[2]^3 * c[4] + c[1]^3 *
                            (c[3] + 10 * c[5]) + 30 * c[2]^2 *
                            (6 * c[3]^2 + 64 * c[4]^2 +
                               140 * c[3] * c[5] + 945 * c[5]^2) +
                            c[1]^2 * (2 * c[2]^2 + 18 * c[3]^2 +
                                        36 * c[2] * c[4] +
                                        192 * c[4]^2 +
                                        375 * c[3] * c[5] +
                                        2250 * c[5]^2) +
                            36 * c[2] * c[4] *
                            (125 * c[3]^2 + 528 * c[4]^2 +
                               3360 * c[3] * c[5] +
                               25725 * c[5]^2) +
                            3 * c[1] *
                            (45 * c[3]^3 + 1584 * c[3] * c[4]^2 +
                               1590 * c[3]^2 * c[5] +
                               21360 * c[4]^2 * c[5] +
                               21525 * c[3] * c[5]^2 +
                               110250 * c[5]^3 +
                               12 * c[2]^2 * (c[3] + 10 * c[5]) +
                               8 * c[2] * c[4] *
                               (32 * c[3] + 375 * c[5])) +
                            9 * (45 * c[3]^4 + 8704 * c[4]^4 +
                                   2415 * c[3]^3 * c[5] +
                                   932400 * c[4]^2 * c[5]^2 +
                                   3018750 * c[5]^4 +
                                   20 * c[3]^2 *
                                   (178 * c[4]^2 +
                                      2765 * c[5]^2) +
                                   35 * c[3] *
                                   (3104 * c[4]^2 * c[5] +
                                      18075 * c[5]^3)))
      return(pskurtosis)
    }
    if (length(xstart) == 0) {
      set.seed(seed)
      cstart1 <- runif(n, min = 0, max = 2)
      cstart2 <- runif(n, min = -1, max = 1)
      cstart3 <- runif(n, min = -1, max = 1)
      cstart4 <- runif(n, min = -0.025, max = 0.025)
      cstart5 <- runif(n, min = -0.025, max = 0.025)
      lstart1 <- runif(n, min = -5, max = 5)
      lstart2 <- runif(n, min = -0.025, max = 0.025)
      lstart3 <- runif(n, min = -0.025, max = 0.025)
      lstart4 <- runif(n, min = -0.025, max = 0.025)
      xstart <- cbind(cstart1, cstart2, cstart3, cstart4,
                      cstart5, lstart1, lstart2, lstart3,
                      lstart4)
    }
    test <- numeric(nrow(xstart))
    for (i in 1:nrow(xstart)) {
      test[i] <-
        ifelse(!is.na(poly_skurt_check(xstart[i, ],
                                       a = c(skews, fifths,
                                             sixths))[1]),
               TRUE, FALSE)
    }
    xstart2 <- xstart[which(test == TRUE), , drop = FALSE]
    constants0 <- NULL
    if (is.null(xstart2)) {
      stop(error1)
    } else {
      converged <- NULL
      for (i in 1:nrow(xstart2)) {
        nl_solution <- nleqslv(x = xstart2[i, ],
                               fn = poly_skurt_check,
                               a = c(skews, fifths, sixths),
                               method = "Broyden",
                               control = list(ftol = 1e-05))
        if (nl_solution$termcd == 1) {
          converged <- rbind(converged, nl_solution$x)
        }
      }
      constants <- converged[!duplicated(converged), , drop = FALSE]
      if (is.null(converged)) {
        constants <- data.frame()
        if (length(Six) == 0) {
          stop(error1)
        } else {
          v <- 1
          while (nrow(constants) == 0) {
            xstart2 <- NULL
            converged <- NULL
            while (is.null(xstart2)) {
              test <- numeric(nrow(xstart))
              for (i in 1:nrow(xstart)) {
                test[i] <-
                  ifelse(!is.na(poly_skurt_check(xstart[i, ],
                                                 a = c(skews,
                                                       fifths,
                                                       sixths +
                                                         Six[v]))[1]),
                         TRUE, FALSE)
              }
              xstart2 <- xstart[which(test == TRUE), , drop = FALSE]
              v <- v + 1
              if (v > length(Six)) break
            }
            if (is.null(xstart2)) {
              stop(error1)
            } else {
              for (i in 1:nrow(xstart2)) {
                nl_solution <- nleqslv(x = xstart2[i, ],
                                       fn = poly_skurt_check,
                                       a = c(skews, fifths,
                                             sixths + Six[v - 1]),
                                       method = "Broyden",
                                       control = list(ftol = 1e-05))
                if (nl_solution$termcd == 1) {
                  converged <- rbind(converged, nl_solution$x)
                }
              }
            }
            constants <- converged[!duplicated(converged), , drop = FALSE]
            if(v > length(Six)) break
          }
          SixCorr1 <- Six[v - 1]
        }
      }
      if (is.null(converged)) {
        stop(error1)
      }
      constants2 <- data.frame()
      if (nrow(constants) != 0) {
        constants <- data.frame(cbind(-constants[, 2] - 3 * constants[, 4],
                                      constants),
                                rep("FALSE", nrow(constants)))
        colnames(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5",
                                 "lambda1", "lambda2", "lambda3",
                                 "lambda4", "valid.pdf")
        constants$valid.pdf <- as.character(constants$valid.pdf)
        for (i in 1:nrow(constants)) {
          if (power_norm_corr(as.numeric(constants[i, 1:6]),
                              method = method) < 0) {
            constants[i, 2] <- -constants[i, 2]
            constants[i, 4] <- -constants[i, 4]
            constants[i, 6] <- -constants[i, 6]
          }
          check <- suppressWarnings(pdf_check(as.numeric(constants[i, 1:6]),
                                            method = method)$valid.pdf)
          if (check[1] == TRUE) {
            constants$valid.pdf[i] <- "TRUE"
          } else {
            constants$valid.pdf[i] <- "FALSE"
          }
        }
        constants$skurtosis <- apply(constants[, 2:6], 1, poly_skurt)
        constants <-  as.data.frame(cbind(rep(skews, nrow(constants)),
                                          rep(fifths, nrow(constants)),
                                          rep(sixths, nrow(constants)),
                                          constants))
        colnames(constants) <- c("skew", "fifth", "sixth", "c0", "c1",
                                 "c2", "c3", "c4", "c5", "lambda1",
                                 "lambda2", "lambda3", "lambda4",
                                 "valid.pdf", "skurtosis")
        constants0 <- constants
        constants2 <- constants[constants$valid.pdf == "TRUE", , drop = FALSE]
      }
      if (nrow(constants2) != 0) {
        m <- min(constants2$skurtosis)
        C2 <- constants2[constants2$skurtosis == m, , drop = FALSE]
        stop.time <- Sys.time()
        Time <- round(difftime(stop.time, start.time,
                               units = "min"), 3)
        return(list(Min = C2, C = constants2, Invalid.C = NULL,
                    start = xstart2, Time = Time,
                    SixCorr1 = SixCorr1, SkurtCorr1 = SkurtCorr1))
      }
      if (nrow(constants2) == 0) {
        if (length(Skurt) == 0) {
          m <- min(constants0$skurtosis)
          C2 <- constants0[constants0$skurtosis == m, , drop = FALSE]
          stop.time <- Sys.time()
          Time <- round(difftime(stop.time, start.time,
                                 units = "min"), 3)
          cat(error2)
          return(list(Min = C2, C = NULL, Invalid.C = constants0,
                      start = xstart2, Time = Time,
                      SixCorr1 = SixCorr1, SkurtCorr1 = SkurtCorr1))
        } else {
          v <- 1
          while (nrow(constants2) == 0) {
            con_solution <-
              suppressWarnings(find_constants(method = method,
                                              skews = skews,
                                              skurts =
                                                min(constants0$skurtosis) +
                                                Skurt[v],
                                              fifths = fifths,
                                              sixths = sixths, Six = NULL,
                                              cstart = NULL,
                                              n = n, seed = seed))
            if (length(con_solution) != 1) {
              if (con_solution$valid == "TRUE") {
                constants2 <-
                  matrix(con_solution$constants, nrow = 1,
                         ncol = length(con_solution$constants))
              }
            }
            v <- v + 1
            if (v > length(Skurt)) break
          }
          SkurtCorr1 <- Skurt[v - 1]
        }
        if (nrow(constants2) != 0) {
          constants2 <-  data.frame(skews, fifths, sixths,
                                    constants2, "TRUE",
                                    min(constants0$skurtosis) +
                                      Skurt[v-1])
          colnames(constants2) <- c("skew", "fifth", "sixth",
                                    "c0", "c1", "c2", "c3", "c4",
                                    "c5", "valid.pdf", "skurtosis")
          m <- min(constants2$skurtosis)
          C2 <- constants2[constants2$skurtosis == m, , drop = FALSE]
          stop.time <- Sys.time()
          Time <- round(difftime(stop.time, start.time,
                                 units = "min"), 3)
          return(list(Min = C2, C = constants2,
                      Invalid.C = constants0, start = xstart2,
                      Time = Time, SixCorr1 = SixCorr1,
                      SkurtCorr1 = SkurtCorr1))
        }
        if (nrow(constants2) == 0) {
          m <- min(constants0$skurtosis)
          C2 <- constants0[constants0$skurtosis == m, , drop = FALSE]
          stop.time <- Sys.time()
          Time <- round(difftime(stop.time, start.time,
                                 units = "min"), 3)
          cat(error2)
          return(list(Min = C2, C = NULL, Invalid.C = constants0,
                      start = xstart2, Time = Time,
                      SixCorr1 = SixCorr1, SkurtCorr1 = SkurtCorr1))
        }
      }
    }
  }
}
