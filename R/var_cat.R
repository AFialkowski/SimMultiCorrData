#' @title Calculate Variance of Binary or Ordinal Variable
#'
#' @description This function calculates the variance of a binary or ordinal (r > 2 categories) variable.  It uses the
#'     formula given by Olsson et al. (1982, \doi{10.1007/BF02294164}) in describing polyserial and point-polyserial correlations.  The
#'     function is used to find intercorrelations involving ordinal variables or variables that are treated as ordinal
#'     (i.e. count variables in the method used in \code{\link[SimMultiCorrData]{rcorrvar2}}).
#'     For an ordinal variable with r >= 2 categories, the variance is given by:
#'     \deqn{\sum_{j=1}^{r} {y_{j}}^{2}*p_{j} - {(\sum_{j=1}^{r} y_{j}*p_{j})}^{2}}.  Here, \eqn{y_{j}} is the j-th support
#'     value and \eqn{p_{j}} is \eqn{Pr(Y = y_{j})}.  This function would not ordinarily be called by the user.
#' @param marginal a vector of cumulative probabilities defining the marginal distribution of the variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param support a vector of containing the ordered support values
#' @export
#' @keywords variance
#' @seealso \code{\link[SimMultiCorrData]{ordnorm}}, \code{\link[SimMultiCorrData]{rcorrvar}},
#'     \code{\link[SimMultiCorrData]{rcorrvar2}}, \code{\link[SimMultiCorrData]{findintercorr_cont_cat}},
#'     \code{\link[SimMultiCorrData]{findintercorr_cont_pois2}}, \cr \code{\link[SimMultiCorrData]{findintercorr_cont_nb2}}
#' @return A scalar equal to the variance
#' @references Olsson U, Drasgow F, & Dorans NJ (1982). The Polyserial Correlation Coefficient. Psychometrika, 47(3): 337-47.
#'     \doi{10.1007/BF02294164}.
#'
var_cat <- function(marginal, support) {
  if (length(marginal) == 1) {
    return(marginal[1] * (1 - marginal[1]))
  } else {
    one <- rep(NA, length(support))
    two <- rep(NA, length(support))
    one[1] <- (support[1]^2) * marginal[1]
    for (i in 2:(length(support) - 1)) {
      one[i] <- (support[i]^2) * (marginal[i] - marginal[i - 1])
    }
    one[length(support)] <- (support[length(support)]^2) *
      (1 - marginal[length(marginal)])
    two[1] <- support[1] * marginal[1]
    for (i in 2:(length(support) - 1)) {
      two[i] <- support[i] * (marginal[i] - marginal[i - 1])
    }
    two[length(support)] <- support[length(support)] *
      (1 - marginal[length(marginal)])
    return(sum(one) - (sum(two)^2))
  }
}
