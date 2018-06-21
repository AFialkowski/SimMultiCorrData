#' @title Calculate Intermediate MVN Correlation for Poisson Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_pois x k_pois} intermediate matrix of correlations for the
#'     Poisson variables using the method of Yahav & Shmueli (2012, \doi{10.1002/asmb.901}). The intermediate correlation between Z1 and Z2 (the
#'     standard normal variables used to generate the Poisson variables Y1 and Y2 via the inverse cdf method) is
#'     calculated using a logarithmic transformation of the target correlation.  First, the upper and lower Frechet-Hoeffding bounds
#'     (mincor, maxcor) \eqn{\rho_{y1,y2}} are simulated.  Then the intermediate correlation is found as follows:
#'     \deqn{\rho_{z1,z2} = (1/b) * log((\rho_{y1,y2} - c)/a)}, where \eqn{a = -(maxcor * mincor)/(maxcor + mincor)},
#'     \eqn{b = log((maxcor + a)/a)}, and \eqn{c = -a}.  The function adapts code from Amatya & Demirtas' (2016) package
#'     \code{\link[PoisNor]{PoisNor-package}} by:
#'
#'     1) allowing specifications for the number of random variates and the seed for reproducibility
#'
#'     2) providing the following checks: if \eqn{\rho_{z1,z2}} >= 1, \eqn{\rho_{z1,z2}} is set to 0.99; if \eqn{\rho_{z1,z2}} <= -1,
#'     \eqn{\rho_{z1,z2}} is set to -0.99.
#'
#'     The function is used in \code{\link[SimMultiCorrData]{findintercorr}} and \code{\link[SimMultiCorrData]{rcorrvar}}.
#'     This function would not ordinarily be called by the user.
#'
#'     Note: The method used here is also used in the packages \code{\link[PoisBinOrdNor]{PoisBinOrdNor-package}} and
#'     \code{\link[PoisBinOrdNonNor]{PoisBinOrdNonNor-package}} by Demirtas et al. (2017), but without my modifications.
#'
#' @param rho_pois a \code{k_pois x k_pois} matrix of target correlations
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{Poisson}})
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @import stats
#' @import utils
#' @export
#' @keywords intermediate, correlation, Poisson, method 1
#' @seealso \code{\link[PoisNor]{PoisNor-package}}, \code{\link[SimMultiCorrData]{findintercorr_nb}},
#'     \code{\link[SimMultiCorrData]{findintercorr_pois_nb}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return the \code{k_pois x k_pois} intermediate correlation matrix for the Poisson variables
#' @references
#' Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15): 3129-39. \doi{10.1080/00949655.2014.953534}.
#'
#' Amatya A & Demirtas H (2016). PoisNor: Simultaneous Generation of Multivariate Data with Poisson and Normal Marginals.
#'     R package version 1.1. \url{https://CRAN.R-project.org/package=PoisNor}
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2): 104-109.
#'
#' Demirtas H, Hu Y, & Allozi R (2017). PoisBinOrdNor: Data Generation with Poisson, Binary, Ordinal and Normal
#'     Components. R package version 1.4. \url{https://CRAN.R-project.org/package=PoisBinOrdNor}
#'
#' Demirtas H, Nordgren R, & Allozi R (2017). PoisBinOrdNonNor: Generation of Up to Four Different
#'     Types of Variables. R package version 1.3. \url{https://CRAN.R-project.org/package=PoisBinOrdNonNor}
#'
#' Frechet M.  Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA.  1951;14:53-77.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1): 91-102. \doi{10.1002/asmb.901}.
#'
findintercorr_pois <- function(rho_pois, lam, nrand = 100000, seed = 1234) {
  set.seed(seed)
  u <- runif(nrand, 0, 1)
  Sigma_pois <- diag(1, nrow(rho_pois), ncol(rho_pois))
  for (i in 1:(nrow(rho_pois) - 1)) {
    for (j in (i + 1):ncol(rho_pois)) {
      maxcor <- cor(qpois(u, lam[i]), qpois(u, lam[j]))
      mincor <- cor(qpois(u, lam[i]), qpois(1 - u, lam[j]))
      a <- -(maxcor * mincor)/(maxcor + mincor)
      b <- log((maxcor + a)/a)
      c <- -a
      Sigma_pois[i, j] <- (1/b) * log((rho_pois[i, j] - c)/a)
      if (Sigma_pois[i, j] >= 1) {
        Sigma_pois[i, j] <- 0.99
      }
      if (Sigma_pois[i, j] <= -1) {
        Sigma_pois[i, j] <- -0.99
      }
      Sigma_pois[j, i] <- Sigma_pois[i, j]
    }
  }
  return(Sigma_pois)
}
