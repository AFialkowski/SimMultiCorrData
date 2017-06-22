#' @title Calculate Upper Frechet-Hoeffding Correlation Bound: Poisson - Normal Variables
#'
#' @description This function calculates the upper Frechet-Hoeffding bound on the correlation between a Poisson variable
#'     and the normal variable used to generate it.  It is used in \code{\link[SimMultiCorrData]{findintercorr_cat_pois}}
#'     and \code{\link[SimMultiCorrData]{findintercorr_cont_pois}} in calculating the intermediate MVN correlations.  This uses
#'     the method of Amatya & Demirtas (2015).  This function would not ordinarily be called directly by the user.
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{dpois}})
#' @param n_unif the number of uniform random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @import stats
#' @import utils
#' @export
#' @keywords correlation, Poisson, method 1
#' @seealso \code{\link[SimMultiCorrData]{findintercorr_cat_pois}}, \code{\link[SimMultiCorrData]{findintercorr_cont_pois}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}
#' @return A scalar equal to the correlation upper bound.
#' @references Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15): 3129-39.
#'
#' Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1): 91-102. \doi{10.1002/asmb.901}.
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2): 104-109.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Frechet M.  Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA.  1951;14:53-77.
chat_pois <- function(lam, n_unif = 10000, seed = 1234) {
  set.seed(seed)
  u <- runif(n_unif, 0, 1)
  chat <- cor(qpois(u, lam), qnorm(u, 0, 1))
  return(chat)
}
