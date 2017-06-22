#' @title Calculate Final Correlation Matrix
#'
#' @description This function calculates the final correlation matrix based on simulated variable type (ordinal, continuous, Poisson, and/or
#'     Negative Binomial).  The function is used in \code{\link[SimMultiCorrData]{rcorrvar}} and
#'     \code{\link[SimMultiCorrData]{rcorrvar2}}.  This would not ordinarily be called directly by the user.
#'
#' @param k_cat the number of ordinal (r >= 2 categories) variables
#' @param k_cont the number of continuous variables
#' @param k_pois the number of Poisson variables
#' @param k_nb the number of Negative Binomial variables
#' @param Y_cat the ordinal (r >= 2 categories) variables
#' @param Yb the continuous variables
#' @param Y_pois the Poisson variables
#' @param Y_nb the Negative Binomial variables
#' @import stats
#' @import utils
#' @export
#' @seealso \code{\link[SimMultiCorrData]{rcorrvar}}, \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @return a correlation matrix
calc_final_corr <- function(k_cat, k_cont, k_pois, k_nb,
                            Y_cat, Yb, Y_pois, Y_nb) {
  if (k_cat > 0 & k_cont == 0 & k_pois == 0 & k_nb == 0) {
    rho_calc <- cor(Y_cat)
  }
  if (k_cat == 0 & k_cont > 0 & k_pois == 0 & k_nb == 0) {
    rho_calc <- cor(Yb)
  }
  if (k_cat == 0 & k_cont == 0 & k_pois > 0 & k_nb == 0) {
    rho_calc <- cor(Y_pois)
  }
  if (k_cat == 0 & k_cont == 0 & k_pois == 0 & k_nb > 0) {
    rho_calc <- cor(Y_nb)
  }
  if (k_cat > 0 & k_cont > 0 & k_pois > 0 & k_nb > 0) {
    rho_calc <- cor(cbind(Y_cat, Yb, Y_pois, Y_nb))
  }
  if (k_cat > 0 & k_cont > 0 & k_pois > 0 & k_nb == 0) {
    rho_calc <- cor(cbind(Y_cat, Yb, Y_pois))
  }
  if (k_cat > 0 & k_cont > 0 & k_pois == 0 & k_nb > 0) {
    rho_calc <- cor(cbind(Y_cat, Yb, Y_nb))
  }
  if (k_cat > 0 & k_cont == 0 & k_pois > 0 & k_nb > 0) {
    rho_calc <- cor(cbind(Y_cat, Y_pois, Y_nb))
  }
  if (k_cat == 0 & k_cont > 0 & k_pois > 0 & k_nb > 0) {
    rho_calc <- cor(cbind(Yb, Y_pois, Y_nb))
  }
  if (k_cat > 0 & k_cont > 0 & k_pois == 0 & k_nb == 0) {
    rho_calc <- cor(cbind(Y_cat, Yb))
  }
  if (k_cat > 0 & k_cont == 0 & k_pois > 0 & k_nb == 0) {
    rho_calc <- cor(cbind(Y_cat, Y_pois))
  }
  if (k_cat > 0 & k_cont == 0 & k_pois == 0 & k_nb > 0) {
    rho_calc <- cor(cbind(Y_cat, Y_nb))
  }
  if (k_cat == 0 & k_cont > 0 & k_pois > 0 & k_nb == 0) {
    rho_calc <- cor(cbind(Yb, Y_pois))
  }
  if (k_cat == 0 & k_cont > 0 & k_pois == 0 & k_nb > 0) {
    rho_calc <- cor(cbind(Yb, Y_nb))
  }
  if (k_cat == 0 & k_cont == 0 & k_pois > 0 & k_nb > 0) {
    rho_calc <- cor(cbind(Y_pois, Y_nb))
  }
  return(rho_calc)
}
