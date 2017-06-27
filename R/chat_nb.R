#' @title Calculate Upper Frechet-Hoeffding Correlation Bound: Negative Binomial - Normal Variables
#'
#' @description This function calculates the upper Frechet-Hoeffding bound on the correlation between a Negative Binomial variable
#'     and the normal variable used to generate it.  It is used in \code{\link[SimMultiCorrData]{findintercorr_cat_nb}}
#'     and \code{\link[SimMultiCorrData]{findintercorr_cont_nb}} in calculating the intermediate MVN correlations.  This extends
#'     the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534}) to Negative Binomial variables.  This function would not ordinarily be called directly by the user.
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param n_unif the number of uniform random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @import stats
#' @import utils
#' @export
#' @keywords correlation, Negative Binomial, method 1
#' @seealso \code{\link[SimMultiCorrData]{findintercorr_cat_nb}}, \code{\link[SimMultiCorrData]{findintercorr_cont_nb}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}
#' @return A scalar equal to the correlation upper bound.
#' @references Please see references for \code{\link[SimMultiCorrData]{chat_pois}}.
#'
chat_nb <- function(size, prob, mu = NULL, n_unif = 10000, seed = 1234) {
  set.seed(seed)
  u <- runif(n_unif, 0, 1)
  if (length(prob) > 0) {
    chat <- cor(qnbinom(u, size, prob), qnorm(u, 0, 1))
  }
  if (length(mu) > 0) {
    chat <- cor(qnbinom(u, size, mu = mu), qnorm(u, 0, 1))
  }
  return(chat)
}
