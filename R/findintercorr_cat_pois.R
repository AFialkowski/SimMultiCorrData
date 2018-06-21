#' @title Calculate Intermediate MVN Correlation for Ordinal - Poisson Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_cat x k_pois} intermediate matrix of correlations for the \code{k_cat} ordinal (r >=
#'     2 categories) and \code{k_pois} Poisson variables. It extends the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534})
#'     to ordinal - Poisson pairs.
#'     Here, the intermediate correlation between Z1 and Z2 (where Z1 is the standard normal variable discretized to produce an
#'     ordinal variable Y1, and Z2 is the standard normal variable used to generate a Poisson variable via the inverse cdf method) is
#'     calculated by dividing the target correlation by a correction factor.  The correction factor is the product of the
#'     upper Frechet-Hoeffding bound on the correlation between a Poisson variable and the normal variable used to generate it
#'     (see \code{\link[SimMultiCorrData]{chat_pois}}) and a simulated GSC upper bound on the correlation between an ordinal variable
#'     and the normal variable used to generate it (see Demirtas & Hedeker, 2011, \doi{10.1198/tast.2011.10090}).  The function is used in
#'     \code{\link[SimMultiCorrData]{findintercorr}} and \code{\link[SimMultiCorrData]{rcorrvar}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param rho_cat_pois a \code{k_cat x k_pois} matrix of target correlations among ordinal and Poisson variables
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{Poisson}})
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @import stats
#' @import utils
#' @export
#' @keywords intermediate, correlation, ordinal, Poisson, method 1
#' @seealso \code{\link[SimMultiCorrData]{chat_pois}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return a k_cat x k_pois matrix whose rows represent the k_cat ordinal variables and columns represent the k_pois Poisson variables
#' @references Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15): 3129-39. \doi{10.1080/00949655.2014.953534}.
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2): 104-109. \doi{10.1198/tast.2011.10090}.
#'
#' Frechet M.  Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA.  1951;14:53-77.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1): 91-102. \doi{10.1002/asmb.901}.
#'
findintercorr_cat_pois <- function(rho_cat_pois, marginal, lam,
                                   nrand = 100000, seed = 1234) {
  Sigma_cat_pois <- matrix(1, nrow = nrow(rho_cat_pois),
                           ncol = ncol(rho_cat_pois))
  set.seed(seed)
  n2 <- rnorm(nrand, 0, 1)
  for (i in 1:nrow(rho_cat_pois)) {
    for (j in 1:ncol(rho_cat_pois)) {
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
      Sigma_cat_pois[i, j] <-
        rho_cat_pois[i, j]/(chat_pois(lam[j], n_unif = nrand, seed = seed) *
                              cor(yord[order(yord)], n2[order(n2)]))
    }
  }
  return(Sigma_cat_pois)
}
