#' @title Calculate Power Method Correlation
#'
#' @description This function calculates the correlation between a continuous variable, Y1, generated using a third or fifth-
#'     order polynomial transformation and the generating standard normal variable, Z1.  The power method correlation
#'     (described in Headrick & Kowalchuk, 2007, \doi{10.1080/10629360600605065}) is given by:
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
#' @references Please see references for \code{\link[SimMultiCorrData]{pdf_check}}.
#'
#' @examples
#' # Beta(a = 4, b = 2) Distribution
#' power_norm_corr(c = c(0.108304, 1.104252, -0.123347, -0.045284, 0.005014,
#'                       0.001285),
#'                 method = "Polynomial")
#'
#' # Switch signs on c1, c3, and c5 to get negative correlation (invalid pdf):
#' power_norm_corr(c = c(0.108304, -1.104252, -0.123347, 0.045284, 0.005014,
#'                       -0.001285),
#'                 method = "Polynomial")
#'
#' @export
power_norm_corr <- function(c, method) {
  c <- as.numeric(c)
  if (method == "Fleishman") return(c[2] + 3 * c[4])
  if (method == "Polynomial") return(c[2] + 3 * c[4] + 15 * c[6])
}

