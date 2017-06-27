#' @title Calculate Denominator Used in Intercorrelations Involving Ordinal Variables
#'
#' @description This function calculates part of the the denominator used to find intercorrelations involving ordinal variables
#'     or variables that are treated as ordinal (i.e. count variables in the method used in
#'     \code{\link[SimMultiCorrData]{rcorrvar2}}).  It uses the formula given by Olsson et al. (1982, \doi{10.1007/BF02294164}) in
#'     describing polyserial and point-polyserial correlations.  For an ordinal variable with r >= 2 categories, the value is given by:
#'     \deqn{\sum_{j = 1}^{r-1} \phi(\tau_{j})*(y_{j+1} - y_{j}),} where
#'     \deqn{\phi(\tau) = (2\pi)^{-1/2} * exp(-0.5 * \tau^2).}  Here, \eqn{y_{j}} is the j-th support
#'     value and \eqn{\tau_{j}} is \eqn{\Phi^{-1}(\sum_{i=1}^{j} Pr(Y = y_{i}))}.  This function would not ordinarily be called directly by the user.
#' @param marginal a vector of cumulative probabilities defining the marginal distribution of the variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param support a vector of containing the ordered support values
#' @import stats
#' @import utils
#' @export
#' @keywords intercorrelation, ordinal
#' @seealso \code{\link[SimMultiCorrData]{ordnorm}}, \code{\link[SimMultiCorrData]{rcorrvar}},
#'     \code{\link[SimMultiCorrData]{rcorrvar2}}, \code{\link[SimMultiCorrData]{findintercorr_cont_cat}},
#'     \code{\link[SimMultiCorrData]{findintercorr_cont_pois2}}, \cr
#'     \code{\link[SimMultiCorrData]{findintercorr_cont_nb2}}
#' @return A scalar
#' @references
#' Olsson U, Drasgow F, & Dorans NJ (1982). The Polyserial Correlation Coefficient. Psychometrika, 47(3): 337-47.
#'     \doi{10.1007/BF02294164}.
#'
denom_corr_cat <- function(marginal, support) {
  denom <- rep(NA, length(marginal))
  for (j in 1:length(marginal)) {
    denom[j] <- (1/sqrt(2 * pi)) * exp(-0.5 * ((qnorm(marginal[j]))^2)) *
      (support[j + 1] - support[j])
  }
  return(sum(denom))
}
