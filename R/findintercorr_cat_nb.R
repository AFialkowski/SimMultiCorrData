#' @title Calculate Intermediate MVN Correlation for Ordinal - Negative Binomial Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_cat x k_nb} intermediate matrix of correlations for the k_cat ordinal (r >=
#'     2 categories) and k_nb Negative Binomial variables. It extends the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534})
#'     to ordinal - Negative Binomial pairs.  Here, the intermediate correlation between Z1 and Z2 (where Z1 is the standard normal variable
#'     discretized to produce an ordinal variable Y1, and Z2 is the standard normal variable used to generate a Negative Binomial
#'     variable via the inverse cdf method) is calculated by dividing the target correlation by a correction factor.  The
#'     correction factor is the product of the upper Frechet-Hoeffding bound on the correlation between a Negative Binomial variable
#'     and the normal variable used to generate it (see \code{\link[SimMultiCorrData]{chat_nb}}) and a simulated GSC upper bound on
#'     the correlation between an ordinal variable and the normal variable used to generate it (see Demirtas & Hedeker, 2011,
#'     \doi{10.1198/tast.2011.10090}).
#'     The function is used in \code{\link[SimMultiCorrData]{findintercorr}} and \code{\link[SimMultiCorrData]{rcorrvar}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param rho_cat_nb a \code{k_cat x k_nb} matrix of target correlations among ordinal and Negative Binomial variables
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{NegBinomial}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @import stats
#' @import utils
#' @export
#' @keywords intermediate, correlation, ordinal, Negative Binomial, method 1
#' @seealso \code{\link[SimMultiCorrData]{chat_nb}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return a \code{k_cat x k_nb} matrix whose rows represent the \code{k_cat} ordinal variables and columns represent the
#'     \code{k_nb} Negative Binomial variables
#' @references Please see references for \code{\link[SimMultiCorrData]{findintercorr_cat_pois}}
#'
findintercorr_cat_nb <- function(rho_cat_nb, marginal, size, prob,
                                 mu = NULL, nrand = 100000, seed = 1234) {
  Sigma_cat_nb <- matrix(1, nrow = nrow(rho_cat_nb), ncol = ncol(rho_cat_nb))
  set.seed(seed)
  n2 <- rnorm(nrand, 0, 1)
  for (i in 1:nrow(rho_cat_nb)) {
    for (j in 1:ncol(rho_cat_nb)) {
      yord <- numeric(length(n2))
      for (r in 1:length(marginal[[i]])) {
        if (r != length(marginal[[i]])) {
          q1 <- qnorm(marginal[[i]][r])
          q2 <- qnorm(marginal[[i]][r + 1])
          yord[(q1 < n2) & (n2 <= q2)] <- r
        } else {
          yord[n2 > qnorm(marginal[[i]][r])] <- r
        }
      }
      yord <- yord + 1
      if (length(prob) > 0) {
        Sigma_cat_nb[i, j] <-
          rho_cat_nb[i, j]/(chat_nb(size[j], prob[j], n_unif = nrand,
                                    seed = seed) * cor(yord[order(yord)],
                                                       n2[order(n2)]))
      }
      if (length(mu) > 0) {
        Sigma_cat_nb[i, j] <-
          rho_cat_nb[i, j]/(chat_nb(size[j], mu = mu[j], n_unif = nrand,
                                    seed = seed) * cor(yord[order(yord)],
                                                       n2[order(n2)]))
      }
    }
  }
  return(Sigma_cat_nb)
}
