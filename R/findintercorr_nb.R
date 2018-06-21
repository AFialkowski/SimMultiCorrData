#' @title Calculate Intermediate MVN Correlation for Negative Binomial Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_nb x k_nb} intermediate matrix of correlations for the Negative Binomial variables by
#'     extending the method of Yahav & Shmueli (2012, \doi{10.1002/asmb.901}). The intermediate correlation between Z1 and Z2 (the
#'     standard normal variables used to generate the Negative Binomial variables Y1 and Y2 via the inverse cdf method) is
#'     calculated using a logarithmic transformation of the target correlation.  First, the upper and lower Frechet-Hoeffding bounds
#'     (mincor, maxcor) on \eqn{\rho_{y1,y2}} are simulated.  Then the intermediate correlation is found as follows:
#'     \deqn{\rho_{z1,z2} = (1/b) * log((\rho_{y1,y2} - c)/a)}, where \eqn{a = -(maxcor * mincor)/(maxcor + mincor)},
#'     \eqn{b = log((maxcor + a)/a)}, and \eqn{c = -a}.  The function adapts code from Amatya & Demirtas' (2016) package
#'     \code{\link[PoisNor]{PoisNor-package}} by:
#'
#'     1) allowing specifications for the number of random variates and the seed for reproducibility
#'
#'     2) providing the following checks: if \eqn{\rho_{z1,z2}} >= 1, \eqn{\rho_{z1,z2}} is set to 0.99; if \eqn{\rho_{z1,z2}} <= -1,
#'     \eqn{\rho_{z1,z2}} is set to -0.99
#'
#'     3) simulating Negative Binomial variables.
#'
#'     The function is used in \code{\link[SimMultiCorrData]{findintercorr}} and \code{\link[SimMultiCorrData]{rcorrvar}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param rho_nb a \code{k_nb x k_nb} matrix of target correlations
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{NegBinomial}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @import stats
#' @import utils
#' @export
#' @keywords intermediate, correlation, Negative Binomial, method 1
#' @seealso \code{\link[PoisNor]{PoisNor-package}}, \code{\link[SimMultiCorrData]{findintercorr_pois}},
#'     \code{\link[SimMultiCorrData]{findintercorr_pois_nb}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return the \code{k_nb x k_nb} intermediate correlation matrix for the Negative Binomial variables
#' @references Please see references for \code{\link[SimMultiCorrData]{findintercorr_pois}}.
#'
findintercorr_nb <- function(rho_nb, size, prob, mu = NULL, nrand = 100000,
                             seed = 1234) {
  set.seed(seed)
  u <- runif(nrand, 0, 1)
  Sigma_nb <- diag(1, nrow(rho_nb), ncol(rho_nb))
  for (i in 1:(nrow(rho_nb) - 1)) {
    for (j in (i + 1):ncol(rho_nb)) {
      if (length(prob) > 0) {
        maxcor <- cor(qnbinom(u, size[i], prob[i]),
                      qnbinom(u, size[j], prob[j]))
        mincor <- cor(qnbinom(u, size[i], prob[i]),
                      qnbinom(1 - u, size[j], prob[j]))
      }
      if (length(mu) > 0) {
        maxcor <- cor(qnbinom(u, size[i], mu = mu[i]),
                      qnbinom(u, size[j], mu = mu[j]))
        mincor <- cor(qnbinom(u, size[i], mu = mu[i]),
                      qnbinom(1 - u, size[j], mu = mu[j]))
      }
      a <- -(maxcor * mincor)/(maxcor + mincor)
      b <- log((maxcor + a)/a)
      c <- -a
      Sigma_nb[i, j] <- (1/b) * log((rho_nb[i, j] - c)/a)
      if (Sigma_nb[i, j] >= 1) {
        Sigma_nb[i, j] <- 0.99
      }
      if (Sigma_nb[i, j] <= -1) {
        Sigma_nb[i, j] <- -0.99
      }
      Sigma_nb[j, i] <- Sigma_nb[i, j]
    }
  }
  return(Sigma_nb)
}
