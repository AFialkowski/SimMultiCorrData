#' @title Find Power Method Transformation Constants
#'
#' @description This function calculates Fleishman's third or Headrick's fifth-order constants necessary to transform a standard normal
#'     random variable into a continuous variable with the specified skewness, standardized kurtosis, and standardized
#'     fifth and sixth cumulants.  It uses \code{\link[BB]{multiStart}} to find solutions to \code{\link[SimMultiCorrData]{fleish}} or
#'     \code{\link[nleqslv]{nleqslv}} for \code{\link[SimMultiCorrData]{poly}}. Multiple starting values are used to ensure the correct
#'     solution is found.  If not user-specified and \code{method} = "Polynomial", the cumulant values are checked to see if they fall in
#'     Headrick's Table 1 (2002, p.691-2, \doi{10.1016/S0167-9473(02)00072-5}) of common distributions (see \code{\link[SimMultiCorrData]{Headrick.dist}}).
#'     If so, his solutions are used as starting values.  Otherwise, a set of \code{n} values randomly generated from uniform distributions is used to
#'     determine the power method constants.
#'
#'     Each set of constants is checked for a positive correlation with the underlying normal variable (using
#'     \code{\link[SimMultiCorrData]{power_norm_corr}})
#'     and a valid power method pdf (using \code{\link[SimMultiCorrData]{pdf_check}}).  If the correlation is <= 0, the signs of c1 and c3 are
#'     reversed (for \code{method} = "Fleishman"), or c1, c3, and c5 (for \code{method} = "Polynomial").  These sign changes have no effect on the cumulants
#'     of the resulting distribution.  If only invalid pdf constants are found and a vector of sixth cumulant correction values (\code{Six}) is provided,
#'     each is checked for valid pdf constants.  The smallest correction that generates a valid power method pdf is used.  If valid pdf constants
#'     still can not be found, the original invalid pdf constants (calculated without a sixth cumulant correction) will be provided if they exist.
#'     If not, the invalid pdf constants calculated with the sixth cumulant correction will be provided.  If no solutions
#'     can be found, an error is given and the result is NULL.
#'
#' @section Reasons for Function Errors:
#'     1) The most likely cause for function errors is that no solutions to \code{\link[SimMultiCorrData]{fleish}} or
#'     \code{\link[SimMultiCorrData]{poly}} converged when using \code{\link[SimMultiCorrData]{find_constants}}.  If this happens,
#'     the simulation will stop.  Possible solutions include: a) increasing the number of initial starting values (\code{n}),
#'     b) using a different seed, or c) specifying a \code{Six} vector of sixth cumulant correction values (for \code{method} = "Polynomial").
#'     If the standardized cumulants are obtained from \code{calc_theory}, the user may need to use rounded values as inputs (i.e.
#'     \code{skews = round(skews, 8)}).  Due to the nature of the integration involved in \code{calc_theory}, the results are
#'     approximations.  Greater accuracy can be achieved by increasing the number of subdivisions (\code{sub}) used in the integration
#'     process.  For example, in order to ensure that skew is exactly 0 for symmetric distributions.
#'
#'     2) In addition, the kurtosis may be outside the region of possible values.  There is an associated lower boundary for kurtosis associated
#'     with a given skew (for Fleishman's method) or skew and fifth and sixth cumulants (for Headrick's method).  Use
#'     \code{\link[SimMultiCorrData]{calc_lower_skurt}} to determine the boundary for a given set of cumulants.
#'
#' @param method the method used to find the constants.  "Fleishman" uses a third-order polynomial transformation and
#'     requires skewness and standardized kurtosis inputs.  "Polynomial" uses Headrick's fifth-order
#'     transformation and requires all four standardized cumulants.
#' @param skews the skewness value
#' @param skurts the standardized kurtosis value (kurtosis - 3, so that normal variables have a value of 0)
#' @param fifths the standardized fifth cumulant (if \code{method} = "Fleishman", keep NULL)
#' @param sixths the standardized sixth cumulant (if \code{method} = "Fleishman", keep NULL)
#' @param Six a vector of correction values to add to the sixth cumulant if no valid pdf constants are found,
#'     ex: \code{Six = seq(1.5, 2,by = 0.05)}; longer vectors take more computation time
#' @param cstart initial value for root-solving algorithm (see \code{\link[BB]{multiStart}} for \code{method} = "Fleishman"
#'     or \code{\link[nleqslv]{nleqslv}} for \code{method} = "Polynomial").  If user-specified,
#'     must be input as a matrix. If NULL and all 4 standardized cumulants (rounded to 3 digits) are within
#'     0.01 of those in Headrick's common distribution table (see \code{\link[SimMultiCorrData]{Headrick.dist}}
#'     data), uses his constants as starting values; else, generates \code{n} sets of random starting values from
#'     uniform distributions.
#' @param n the number of initial starting values to use with root-solver.  More starting values
#'     require more calculation time (default = 25).
#' @param seed the seed value for random starting value generation (default = 1234)
#' @export
#' @import stats
#' @import utils
#' @import BB
#' @import nleqslv
#' @keywords constants, Fleishman, Headrick
#' @seealso \code{\link[BB]{multiStart}}, \code{\link[nleqslv]{nleqslv}},
#'     \code{\link[SimMultiCorrData]{fleish}}, \code{\link[SimMultiCorrData]{poly}},
#'     \code{\link[SimMultiCorrData]{power_norm_corr}}, \code{\link[SimMultiCorrData]{pdf_check}}
#' @return A list with components:
#' @return \code{constants} a vector of valid or invalid power method solutions, c("c0","c1","c2","c3") for \code{method} = "Fleishman" or
#'     c("c0","c1","c2","c3","c4,"c5") for \code{method} = "Polynomial"
#' @return \code{valid} "TRUE" if the constants produce a valid power method pdf, else "FALSE"
#' @return \code{SixCorr1} if \code{Six} is specified, the sixth cumulant correction required to achieve a valid pdf
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
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' Varadhan R, Gilbert PD (2009). BB: An R Package for Solving a Large System of Nonlinear Equations and for
#'     Optimizing a High-Dimensional Nonlinear Objective Function, J. Statistical Software, 32:4,
#'     \url{http://www.jstatsoft.org/v32/i04/}
#'
#' @examples \dontrun{
#' # Compute third-order power method constants.
#'
#' options(scipen = 999) # turn off scientific notation
#'
#' # Exponential Distribution
#' find_constants("Fleishman", 2, 6)
#'
#' # Laplace Distribution
#' find_constants("Fleishman", 0, 3)
#'
#' # Compute fifth-order power method constants.
#'
#' # Logistic Distribution
#' find_constants(method = "Polynomial", skews = 0, skurts = 6/5, fifths = 0,
#'                sixths = 48/7)
#'
#' # with correction to sixth cumulant
#' find_constants(method = "Polynomial", skews = 0, skurts = 6/5, fifths = 0,
#'                sixths = 48/7, Six = seq(1.7, 2, by = 0.01))
#'
#' }
find_constants <- function(method = c("Fleishman", "Polynomial"), skews = NULL,
                           skurts = NULL, fifths = NULL, sixths = NULL,
                           Six = NULL, cstart = NULL, n = 25, seed = 1234) {
  if (method == "Fleishman") {
    error <- "No valid power method constants could be found for the specified
             cumulants.  Try using a different seed or increasing n.\n"
    if (skews == 0 & skurts == 0) {
      constants <- c(0, 1, 0, 0)
      names(constants) <- c("c0", "c1", "c2", "c3")
      return(list(constants = constants, valid = "TRUE"))
      } else {
      if (length(cstart) == 0) {
        set.seed(seed)
        cstart1 <- runif(n, min = -2, max = 2)
        cstart2 <- runif(n, min = -1, max = 1)
        cstart3 <- runif(n, min = -0.5, max = 0.5)
        cstart <- cbind(cstart1, cstart2, cstart3)
      }
      zeros <- multiStart(par = cstart, fn = fleish, gr = NULL,
                          action = "solve", method = c(2, 3, 1),
                          lower = -Inf, upper = Inf, project = NULL,
                          projectArgs = NULL,
                          control = list(trace = FALSE, tol = 1e-05),
                          quiet = TRUE, a = c(skews, skurts))
      converged <- matrix(zeros$par[zeros$converged == "TRUE", ],
                          nrow = sum(zeros$converged == "TRUE"))
      if (nrow(converged) == 0) {
        warning(error)
      } else {
        constants <- converged[!duplicated(converged), , drop = FALSE]
        constants <- data.frame(cbind(-constants[, 2], constants),
                                rep("FALSE", nrow(constants)))
        colnames(constants) <- c("c0", "c1", "c2", "c3", "valid")
        constants$valid <- as.character(constants$valid)
        for (i in 1:nrow(constants)) {
          if (power_norm_corr(as.numeric(constants[i, 1:4]),
                              method = method) < 0) {
            constants[i, 2] <- -constants[i, 2]
            constants[i, 4] <- -constants[i, 4]
          }
          check <- suppressWarnings(pdf_check(as.numeric(constants[i, 1:4]),
                                            method = method)$valid.pdf)
          if (check[1] == TRUE) {
            constants$valid[i] <- "TRUE"
          } else {
            constants$valid[i] <- "FALSE"
          }
        }
        constants2 <- constants[constants$valid == "TRUE", , drop = FALSE]
        if (nrow(constants2) == 0) {
          freq <- table(round(constants[, 2], 6))
          constants <-
            subset(constants, round(constants[, 2], 6) ==
                     as.numeric(names(freq[freq == max(freq)])[1]))
          constants <- as.numeric(constants[1, 1:4])
          names(constants) <- c("c0", "c1", "c2", "c3")
          return(list(constants = constants, valid = "FALSE"))
        } else {
          freq <- table(round(constants2[, 2], 6))
          constants2 <-
            subset(constants2, round(constants2[, 2], 6) ==
                     as.numeric(names(freq[freq == max(freq)])[1]))
          constants2 <- as.numeric(constants2[1, 1:4])
          names(constants2) <- c("c0", "c1", "c2", "c3")
          return(list(constants = constants2, valid = "TRUE"))
        }
      }
    }
  }
  if (method == "Polynomial") {
    SixCorr1 <- NULL
    error <- "No valid power method constants could be found for the specified
             cumulants.  Try using more starting values or increasing the
             value of the sixth cumulant by specifying a Six vector of
             correction values.\n"
    if (skews == 0 & skurts == 0 & fifths == 0 & sixths == 0) {
      constants <- c(0, 1, 0, 0, 0, 0)
      names(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5")
      return(list(constants = constants, valid = "TRUE", SixCorr1 = SixCorr1))
    } else {
      if (length(cstart) == 0) {
        H_dist <- SimMultiCorrData::Headrick.dist
        H <- H_dist[-c(1:4), ]
        M <- H_dist[1:4, ]
        for (m in 1:ncol(M)) {
          if ((abs(round(skews, 3) - round(M[1, m], 3)) < 0.01) &
              (abs(round(skurts, 3) - round(M[2, m], 3)) < 0.01) &
              (abs(round(fifths, 3) - round(M[3, m], 3)) < 0.01) &
              (abs(round(sixths, 3) - round(M[4, m], 3)) < 0.01)) {
            cstart <- matrix(H[2:6, m], nrow = 1, ncol = 5)
          } else {
            cstart <- cstart
          }
        }
        if (length(cstart) == 0) {
          set.seed(seed)
          cstart1 <- runif(n, min = -2, max = 2)
          cstart2 <- runif(n, min = -1, max = 1)
          cstart3 <- runif(n, min = -1, max = 1)
          cstart4 <- runif(n, min = -0.025, max = 0.025)
          cstart5 <- runif(n, min = -0.025, max = 0.025)
          cstart <- cbind(cstart1, cstart2, cstart3, cstart4, cstart5)
        }
      }
      converged <- NULL
      for (i in 1:nrow(cstart)) {
        nl_solution <- nleqslv(x = cstart[i, ], fn = poly,
                               a = c(skews, skurts, fifths, sixths),
                               method = "Broyden",
                               control = list(ftol = 1e-05))
        if (nl_solution$termcd == 1) {
          converged <- rbind(converged, nl_solution$x)
        }
      }
      constants2 <- data.frame()
      constants <- converged[!duplicated(converged), , drop = FALSE]
      constants0 <- NULL
      if (!is.null(converged)) {
        constants <- data.frame(cbind(-constants[, 2] - 3 * constants[, 4],
                                      constants),
                                rep("FALSE", nrow(constants)))
        colnames(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5", "valid")
        constants$valid <- as.character(constants$valid)
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
            constants$valid[i] <- "TRUE"
          } else {
            constants$valid[i] <- "FALSE"
          }
        }
        constants2 <- constants[constants$valid == "TRUE", , drop = FALSE]
        freq <- table(round(constants[, 2], 6))
        constants <- subset(constants, round(constants[, 2], 6) ==
                              as.numeric(names(freq[freq == max(freq)])[1]))
        constants <- as.numeric(constants[1, 1:6])
        names(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5")
        constants0 <- constants
      }
      if (nrow(constants2) == 0) {
        if (length(Six) != 0) {
          v <- 1
          converged <- NULL
          while (nrow(constants2) == 0) {
            for (i in 1:nrow(cstart)) {
              nl_solution <- nleqslv(x = cstart[i, ], fn = poly,
                                     a = c(skews, skurts, fifths,
                                           sixths + Six[v]),
                                     method = "Broyden",
                                     control = list(ftol = 1e-05))
              if (nl_solution$termcd == 1) {
                converged <- rbind(converged, nl_solution$x)
              }
            }
            if (!is.null(converged)) {
              constants <- converged[!duplicated(converged), , drop = FALSE]
              constants <- cbind(-constants[, 2] - 3 * constants[, 4],
                                 constants)
              for (i in 1:nrow(constants)) {
                if (power_norm_corr(as.numeric(constants[i, ]),
                                    method = method) < 0) {
                  constants[i, 2] <- -constants[i, 2]
                  constants[i, 4] <- -constants[i, 4]
                  constants[i, 6] <- -constants[i, 6]
                }
                check <-
                  suppressWarnings(pdf_check(as.numeric(constants[i, ]),
                                           method = method)$valid.pdf)
                if (check[1] == TRUE) {
                  constants2 <- rbind(constants2, constants[i, ])
                } else {
                  constants2 <- constants2
                }
              }
            }
            v <- v + 1
            if (v > length(Six)) break
          }
          SixCorr1 <- Six[v - 1]
        }
      }
      if (nrow(constants2) == 0) {
        if (!is.null(constants0)) {
          return(list(constants = constants0, valid = "FALSE",
                      SixCorr1 = SixCorr1))
        } else {
          if (!is.null(constants)) {
            freq <- table(round(constants[, 2], 6))
            constants <-
              subset(constants, round(constants[, 2], 6) ==
                       as.numeric(names(freq[freq == max(freq)])[1]))
            constants <- as.numeric(constants[1, 1:6])
            names(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5")
            return(list(constants = constants, valid = "FALSE",
                        SixCorr1 = SixCorr1))
          } else {
            warning(error)
          }
        }
      } else {
        constants2 <- as.data.frame(constants2[, 1:6])
        freq <- table(round(constants2[, 2], 6))
        constants2 <-
          subset(constants2, round(constants2[, 2], 6) ==
                   as.numeric(names(freq[freq == max(freq)])[1]))[1, ]
        constants2 <- as.numeric(constants2[1, 1:6])
        names(constants2) <- c("c0", "c1", "c2", "c3", "c4", "c5")
        return(list(constants = constants2, valid = "TRUE",
                    SixCorr1 = SixCorr1))
      }
    }
  }
}
