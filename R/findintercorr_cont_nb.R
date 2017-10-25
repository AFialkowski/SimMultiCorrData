#' @title Calculate Intermediate MVN Correlation for Continuous - Negative Binomial Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_cont x k_nb} intermediate matrix of correlations for the \code{k_cont} continuous and
#'     \code{k_nb} Negative Binomial variables. It extends the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534}) to
#'     continuous variables generated using
#'     Headrick's fifth-order polynomial transformation and Negative Binomial variables.  Here, the intermediate correlation
#'     between Z1 and Z2 (where Z1 is the standard normal variable transformed using Headrick's fifth-order or Fleishman's
#'     third-order method to produce a continuous variable Y1, and Z2 is the standard normal variable used to generate a
#'     Negative Binomial variable via the inverse cdf method) is calculated by dividing the target correlation by a correction factor.
#'     The correction factor is the product of the upper Frechet-Hoeffding bound on the correlation between a Negative Binomial
#'     variable and the normal variable used to generate it (see \code{\link[SimMultiCorrData]{chat_nb}}) and the power method
#'     correlation (described in Headrick & Kowalchuk, 2007, \doi{10.1080/10629360600605065}) between Y1 and Z1.  The function is used in
#'     \code{\link[SimMultiCorrData]{findintercorr}} and \code{\link[SimMultiCorrData]{rcorrvar}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param method the method used to generate the \code{k_cont} continuous variables.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param constants a matrix with \code{k_cont} rows, each a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or
#'     c0, c1, c2, c3, c4, c5 (if \code{method} = "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param rho_cont_nb a \code{k_cont x k_nb} matrix of target correlations among continuous and Negative Binomial variables
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @export
#' @keywords intermediate, correlation, continuous, Fleishman, Headrick, Negative Binomial, method 1
#' @seealso \code{\link[SimMultiCorrData]{chat_nb}}, \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{find_constants}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return a \code{k_cont x k_nb} matrix whose rows represent the \code{k_cont} continuous variables and columns represent the
#'     \code{k_nb} Negative Binomial variables
#' @references Please see references for \code{\link[SimMultiCorrData]{findintercorr_cont_pois}}.
#'
findintercorr_cont_nb <- function(method, constants, rho_cont_nb, size, prob,
                                  mu = NULL, nrand = 100000, seed = 1234) {
  Sigma_cont_nb <- matrix(1, nrow = nrow(rho_cont_nb),
                          ncol = ncol(rho_cont_nb))
  if (length(prob) > 0) {
    for (i in 1:nrow(rho_cont_nb)) {
      for (j in 1:ncol(rho_cont_nb)) {
        Sigma_cont_nb[i, j] <-
          rho_cont_nb[i, j]/(chat_nb(size[j], prob[j], n_unif = nrand,
                                     seed = seed) *
                               power_norm_corr(constants[i, ],
                                               method))
      }
    }
  } else {
    for (i in 1:nrow(rho_cont_nb)) {
      for (j in 1:ncol(rho_cont_nb)) {
        Sigma_cont_nb[i, j] <-
          rho_cont_nb[i, j]/(chat_nb(size[j], mu = mu[j], n_unif = nrand,
                                     seed = seed) *
                               power_norm_corr(constants[i, ],
                                               method))
      }
    }
  }
  return(Sigma_cont_nb)
}
