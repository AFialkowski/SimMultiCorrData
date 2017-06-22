#' @title Separate Target Correlation Matrix by Variable Type
#'
#' @description This function separates the target correlation matrix rho by variable type (ordinal, continuous, Poisson, and/or
#'     Negative Binomial).  The function is used in \code{\link[SimMultiCorrData]{findintercorr}},
#'     \code{\link[SimMultiCorrData]{rcorrvar}}, \code{\link[SimMultiCorrData]{findintercorr2}}, and
#'     \code{\link[SimMultiCorrData]{rcorrvar2}}.  This would not ordinarily be called directly by the user.
#'
#' @param k_cat the number of ordinal (r >= 2 categories) variables
#' @param k_cont the number of continuous variables
#' @param k_pois the number of Poisson variables
#' @param k_nb the number of Negative Binomial variables
#' @param rho the target correlation matrix
#' @export
#' @seealso \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}},
#'     \code{\link[SimMultiCorrData]{findintercorr2}}, \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @return a list containing the target correlation matrix components by variable combination
separate_rho <- function(k_cat, k_cont, k_pois, k_nb, rho) {
  rho_cat <- NA
  rho_cat_cont <- NA
  rho_cat_pois <- NA
  rho_cat_nb <- NA
  rho_cont_cat <- NA
  rho_cont <- NA
  rho_cont_pois <- NA
  rho_cont_nb <- NA
  rho_pois_cat <- NA
  rho_pois_cont <- NA
  rho_pois <- NA
  rho_pois_nb <- NA
  rho_nb_cat <- NA
  rho_nb_cont <- NA
  rho_nb_pois <- NA
  rho_nb <- NA
  if (k_cat > 0) {
    rho_cat <- matrix(rho[1:k_cat, 1:k_cat], nrow = k_cat, ncol = k_cat,
                      byrow = F)
  }
  if (k_cont > 0) {
    rho_cont <- matrix(rho[(k_cat + 1):(k_cat + k_cont),
                           (k_cat + 1):(k_cat + k_cont)],
                       nrow = k_cont, ncol = k_cont, byrow = F)
  }
  if (k_pois > 0) {
    rho_pois <- matrix(rho[(k_cat + k_cont + 1):(k_cat + k_cont + k_pois),
                           (k_cat + k_cont + 1):(k_cat + k_cont + k_pois)],
                       nrow = k_pois, ncol = k_pois, byrow = F)
  }
  if (k_nb > 0) {
    rho_nb <- matrix(rho[(k_cat + k_cont + k_pois + 1):(k_cat + k_cont +
                                                          k_pois + k_nb),
                         (k_cat + k_cont + k_pois + 1):(k_cat + k_cont +
                                                          k_pois + k_nb)],
                     nrow = k_nb, ncol = k_nb, byrow = F)
  }
  if (k_cat > 0 & k_cont > 0) {
    rho_cont_cat <- matrix(rho[(k_cat + 1):(k_cat + k_cont), 1:k_cat],
                           nrow = k_cont, ncol = k_cat, byrow = F)
    rho_cat_cont <- t(rho_cont_cat)
  }
  if (k_cat > 0 & k_pois > 0) {
    rho_cat_pois <- matrix(rho[1:k_cat, (k_cat + k_cont + 1):(k_cat + k_cont +
                                                                k_pois)],
                           nrow = k_cat, ncol = k_pois, byrow = F)
    rho_pois_cat <- t(rho_cat_pois)
  }
  if (k_cat > 0 & k_nb > 0) {
    rho_cat_nb <- matrix(rho[1:k_cat,
                             (k_cat + k_cont + k_pois + 1):(k_cat + k_cont +
                                                              k_pois + k_nb)],
                         nrow = k_cat, ncol = k_nb, byrow = F)
    rho_nb_cat <- t(rho_cat_nb)
  }
  if (k_cont > 0 & k_pois > 0) {
    rho_cont_pois <- matrix(rho[(k_cat + 1):(k_cat + k_cont),
                                (k_cat + k_cont + 1):(k_cat + k_cont +
                                                        k_pois)],
                            nrow = k_cont, ncol = k_pois, byrow = F)
    rho_pois_cont <- t(rho_cont_pois)
  }
  if (k_cont > 0 & k_nb > 0) {
    rho_cont_nb <- matrix(rho[(k_cat + 1):(k_cat + k_cont),
                              (k_cat + k_cont + k_pois + 1):(k_cat +
                                                               k_cont +
                                                               k_pois + k_nb)],
                          nrow = k_cont, ncol = k_nb, byrow = F)
    rho_nb_cont <- t(rho_cont_nb)
  }
  if (k_pois > 0 & k_nb > 0) {
    rho_pois_nb <- matrix(rho[(k_cat + k_cont + 1):(k_cat + k_cont + k_pois),
                              (k_cat + k_cont + k_pois + 1):(k_cat + k_cont +
                                                               k_pois + k_nb)],
                          nrow = k_pois, ncol = k_nb, byrow = F)
    rho_nb_pois <- t(rho_pois_nb)
  }
  return(list(rho_cat = rho_cat, rho_cat_cont = rho_cat_cont,
              rho_cat_pois = rho_cat_pois, rho_cat_nb = rho_cat_nb,
              rho_cont_cat = rho_cont_cat, rho_cont = rho_cont,
              rho_cont_pois = rho_cont_pois, rho_cont_nb = rho_cont_nb,
              rho_pois_cat = rho_pois_cat, rho_pois_cont = rho_pois_cont,
              rho_pois = rho_pois, rho_pois_nb = rho_pois_nb,
              rho_nb_cat = rho_nb_cat, rho_nb_cont = rho_nb_cont,
              rho_nb_pois = rho_nb_pois, rho_nb = rho_nb))
}
