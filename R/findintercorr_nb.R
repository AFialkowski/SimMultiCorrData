#' @title Calculate Intermediate MVN Correlation for Negative Binomial Variables: Method 1
#'
#' @description This function calculates a \code{k_nb x k_nb} intermediate matrix of correlations for the Negative Binomial variables by
#'     extending the method of Yahav & Shmueli (2012). The intermediate correlation between Z1 and Z2 (the
#'     standard normal variables used to generate the Negative Binomial variables Y1 and Y2 via the inverse cdf method) is
#'     calculated using a logarithmic transformation of the target correlation.  First, the upper and lower Frechet-Hoeffding bounds
#'     (mincor, maxcor) on \eqn{\rho_{y1,y2}} are simulated.  Then the intermediate correlation is found as follows:
#'     \deqn{\rho_{z1,z2} = (1/b) * log((\rho_{y1,y2} - c)/a)}, where \eqn{a = -(maxcor * mincor)/(maxcor + mincor)},
#'     \eqn{b = log((maxcor + a)/a)}, and \eqn{c = -a}.  The function adapts code from Amatya & Demirtas' (2016) package
#'     \code{\link[PoisNor]{PoisNor}} by:
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
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @import stats
#' @import utils
#' @export
#' @keywords intermediate, correlation, Negative Binomial, method 1
#' @seealso \code{\link[PoisNor]{PoisNor}}, \code{\link[SimMultiCorrData]{findintercorr_pois}},
#'     \code{\link[SimMultiCorrData]{findintercorr_pois_nb}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return the \code{k_nb x k_nb} intermediate correlation matrix for the Negative Binomial variables
#' @references Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1): 91-102. \doi{10.1002/asmb.901}.
#'
#' Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15): 3129-39.
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2): 104-109.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Frechet M.  Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA.  1951;14:53-77.
#'
#' Amatya A & Demirtas H (2016). PoisNor: Simultaneous Generation of Multivariate Data with Poisson and Normal Marginals.
#'     R package version 1.1. \url{https://CRAN.R-project.org/package=PoisNor}
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
