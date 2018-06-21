#' @title Calculate Maximum Support Value for Count Variables: Correlation Method 2
#'
#' @description This function calculates the maximum support value for count variables by extending the method of Barbiero &
#'     Ferrari (2015, \doi{10.1002/asmb.2072}) to include Negative Binomial variables.  In order for count variables to be treated as ordinal in the
#'     calculation of the intermediate MVN correlation matrix, their infinite support must be truncated (made finite).  This is
#'     done by setting the total cumulative probability equal to 1 - a small user-specified value (\code{pois_eps} or \code{nb_eps}.  The
#'     maximum support value equals the inverse cdf applied to this result.  The values pois_eps and nb_eps may differ for each variable.
#'     The function is used in \code{\link[SimMultiCorrData]{findintercorr2}} and \code{\link[SimMultiCorrData]{rcorrvar2}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param k_pois the number of Poisson variables
#' @param k_nb the number of Negative Binomial variables
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{Poisson}})
#' @param pois_eps a vector of length k_pois containing the truncation values (i.e. = rep(0.0001, k_pois); default = NULL)
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{NegBinomial}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param nb_eps a vector of length \code{k_nb} containing the truncation values (i.e. = rep(0.0001, k_nb); default = NULL)
#' @import stats
#' @import utils
#' @export
#' @keywords intermediate, correlation, Poisson, Negative Binomial, method 2
#' @seealso \code{\link[SimMultiCorrData]{findintercorr2}}, \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @return a data.frame with \code{k_pois + k_nb} rows; the column names are:
#' @return \code{Distribution} Poisson or Negative Binomial
#' @return \code{Number} the variable index
#' @return \code{Max} the maximum support value
#' @references
#' Barbiero A & Ferrari PA (2015). Simulation of correlated Poisson variables. Applied Stochastic Models in
#'     Business and Industry, 31: 669-80. \doi{10.1002/asmb.2072}.
#'
#' Ferrari PA, Barbiero A (2012). Simulating ordinal data, Multivariate Behavioral Research, 47(4): 566-589. \doi{10.1080/00273171.2012.692630}.
#'
max_count_support <- function(k_pois, k_nb, lam, pois_eps = NULL,
                              size, prob, mu = NULL, nb_eps = NULL) {
  max_support <- matrix(1, nrow = k_pois + k_nb, ncol = 2)
  if (k_pois > 0) {
    for (i in 1:k_pois) {
      max_support[i, ] <- append(i, qpois(1 - pois_eps[i], lam[i]))
    }
  }
  if (k_nb > 0) {
    for (i in (k_pois + 1):(k_pois + k_nb)) {
      if (length(prob) > 0) {
        max_support[i, ] <- append(i, qnbinom(1 - nb_eps[i - k_pois],
                                              size[i - k_pois],
                                              prob[i - k_pois]))
      }
      if (length(mu) > 0) {
        max_support[i, ] <- append(i, qnbinom(1 - nb_eps[i - k_pois],
                                              size[i - k_pois],
                                              mu = mu[i - k_pois]))
      }
    }
  }
  max_support <- cbind(append(rep("Poisson", k_pois),
                              rep("Neg_Bin", k_nb)),
                       as.data.frame(max_support))
  colnames(max_support) <- c("Distribution", "Number", "Max")
  return(max_support)
}
