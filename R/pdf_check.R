#' @title Check Polynomial Transformation Constants for Valid Power Method PDF
#'
#' @description This function determines if a given set of constants, calculated using Fleishman's third-order or
#'     Headrick's fifth-order transformation, yields a valid pdf.  This requires 1) the correlation between the
#'     resulting continuous variable and the underlying standard normal variable (see
#'     \code{\link[SimMultiCorrData]{power_norm_corr}}) is > 0, and 2) the constants satisfy certain constraints.
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param method the method used to find the constants.  "Fleishman" uses a third-order polynomial transformation and
#'     "Polynomial" uses Headrick's fifth-order transformation.
#' @export
#' @keywords constants, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{fleish}}, \code{\link[SimMultiCorrData]{poly}},
#'     \code{\link[SimMultiCorrData]{power_norm_corr}}, \code{\link[SimMultiCorrData]{find_constants}}
#' @return A list with components:
#' @return \code{rho_pZ} the correlation between the continuous variable and the underlying standard normal variable
#' @return \code{valid.pdf} "TRUE" if the constants produce a valid power method pdf, else "FALSE"
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
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35.
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3, 65-71.
#' @examples \dontrun{
#' # Chi-squared (df = 1) Distribution (invalid power method pdf)
#' con <- find_constants(method = "Polynomial", skews = sqrt(8), skurts = 12,
#'                       fifths = 48*sqrt(2), sixths = 480)$constants
#' pdf_check(c = con, method = "Polynomial")
#'
#' # Beta (a = 4, b = 2) Distribution (valid power method pdf)
#' con <- find_constants(method = "Polynomial", skews = -0.467707,
#'                       skurts = -0.375, fifths = 1.403122,
#'                       sixths = -0.426136)$constants
#' pdf_check(c = con, method = "Polynomial")
#'
#' }
pdf_check <- function(c, method) {
  c <- as.numeric(c)
  if (method == "Fleishman") {
    c0 <- c[1]
    c1 <- c[2]
    c2 <- c[3]
    c3 <- c[4]
    rho_pZ <- c1 + 3 * c3
    if (rho_pZ == 1) return(list(rho_pZ = rho_pZ, valid.pdf = TRUE))
    if (rho_pZ <= 0) {
      warning("Error: The correlation between Y and Z is non-positive.
          Reverse the signs of c1 and c3.\n")
      valid.pdf <- FALSE
    }
    if (rho_pZ > 1) {
      warning("Error: The correlation between Y and Z is greater than 1.\n")
      valid.pdf <- FALSE
    }
    if ((c3 > (sqrt(5 + 7 * c1^2)/(5 * sqrt(3)) - (2 * c1)/5)) &
        c1 > 0 & c1 < 1 & rho_pZ > 0 & rho_pZ <= 1) {
      warning("This is a valid pdf.\n")
      valid.pdf <- TRUE
    } else {
      warning("This is NOT a valid pdf.\n")
      valid.pdf <- FALSE
    }
    return(list(rho_pZ = rho_pZ, valid.pdf = valid.pdf))
  }
  if (method == "Polynomial") {
    c0 <- c[1]
    c1 <- c[2]
    c2 <- c[3]
    c3 <- c[4]
    c4 <- c[5]
    c5 <- c[6]
    rho_pZ <- c1 + 3 * c3 + 15 * c5
    if (rho_pZ == 1) return(list(rho_pZ = rho_pZ, valid.pdf = TRUE))
    if (rho_pZ <= 0) {
      warning("Error: The correlation between Y and Z is non-positive.
          Reverse the signs of c1, c3, c5.\n")
      valid.pdf <- FALSE
    }
    if (rho_pZ > 1) {
      warning("Error: The correlation between Y and Z is greater than 1.\n")
      valid.pdf <- FALSE
    }
    if (round(c0, 6) == 0 & round(c2, 6) == 0 & round(c4, 6) == 0) {
      z1 <-  sqrt(as.complex((sqrt(as.complex(9 * c3^2 - 20 * c1 * c5)) -
                                3 * c3)/(10 * c5)))
      z2 <-  sqrt(as.complex((-sqrt(as.complex(9 * c3^2 - 20 * c1 * c5)) -
                                3 * c3)/(10 * c5)))
      z3 <- -sqrt(as.complex((sqrt(as.complex(9 * c3^2 - 20 * c1 * c5)) -
                                3 * c3)/(10 * c5)))
      z4 <- -sqrt(as.complex((-sqrt(as.complex(9 * c3^2 - 20 * c1 * c5)) -
                                3 * c3)/(10 * c5)))
    } else {
      S1 <- 54 * c3^3 - 216 * c2 * c3 * c4 + 432 * c1 * c4^2 +
        540 * c5 * c2^2 - 1080 * c1 * c3 * c5
      S2 <- S1^2 - 4 * (9 * c3^2 - 24 * c2 * c4 + 60 * c1 * c5)^3
      S3 <- 3 * c3^2 - 8 * c2 * c4 + 20 * c1 * c5
      S4 <- (4 * c4^2)/(25 * c5^2) + ((sqrt(as.complex(S2)) + S1)^(1/3))/(15 *
             c5 * 2^(1/3)) + (S3 * 2^(1/3))/(5 * c5 * (sqrt(as.complex(S2)) +
             S1)^(1/3)) - (2 * c3)/(5 * c5)
      S5 <- (1/sqrt(as.complex(S4))) * ((48 * c3 * c4)/(100 * c5^2) -
             (16 * c2)/(20 * c5) - (64 * c4^3)/(500 * c5^3))
      S6 <- (8 * c4^2)/(25 * c5^2) - ((sqrt(as.complex(S2)) + S1)^(1/3))/(15 *
             c5 * 2^(1/3)) - (S3 * 2^(1/3))/(5 * c5 * (sqrt(as.complex(S2)) +
             S1)^(1/3)) - (4 * c3)/(5 * c5)
      z1 <- 0.5 * sqrt(as.complex(S4)) + 0.5 * sqrt(as.complex(S6 + S5)) -
        c4/(5 * c5)
      z2 <- -0.5 * sqrt(as.complex(S4)) - 0.5 * sqrt(as.complex(S6 - S5)) -
        c4/(5 * c5)
      z3 <- 0.5 * sqrt(as.complex(S4)) - 0.5 * sqrt(as.complex(S6 + S5)) -
        c4/(5 * c5)
      z4 <- -0.5 * sqrt(as.complex(S4)) + 0.5 * sqrt(as.complex(S6 - S5)) -
        c4/(5 * c5)
    }
    if (rho_pZ > 1 | rho_pZ <= 0 | abs(Im(z1)) < 1e-06 | abs(Im(z2)) < 1e-06 |
        abs(Im(z3)) < 1e-06 | abs(Im(z4)) < 1e-06) {
      warning("This is NOT a valid pdf.  A possible solution is to increase the
           value of the sixth standardized cumulant.\n")
      valid.pdf <- FALSE
    } else {
      warning("This is a valid pdf.\n")
      valid.pdf <- TRUE
    }
    return(list(rho_pZ = rho_pZ, Zs = data.frame(z1 = z1, z2 = z2, z3 = z3,
                                                 z4 = z4),
                valid.pdf = valid.pdf))
  }
}
