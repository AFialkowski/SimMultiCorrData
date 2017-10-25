#' @title Calculate Intermediate MVN Correlation for Continuous - Poisson Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_cont x k_pois} intermediate matrix of correlations for the \code{k_cont} continuous and
#'     \code{k_pois} Poisson variables. It extends the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534}) to continuous
#'     variables generated using Headrick's fifth-order polynomial transformation.  Here, the intermediate correlation between Z1 and Z2 (where Z1 is
#'     the standard normal variable transformed using Headrick's fifth-order or Fleishman's third-order method to produce a
#'     continuous variable Y1, and Z2 is the standard normal variable used to generate a Poisson variable via the inverse cdf method) is
#'     calculated by dividing the target correlation by a correction factor.  The correction factor is the product of the
#'     upper Frechet-Hoeffding bound on the correlation between a Poisson variable and the normal variable used to generate it
#'     (see \code{\link[SimMultiCorrData]{chat_pois}}) and the power method correlation (described in Headrick & Kowalchuk, 2007,
#'     \doi{10.1080/10629360600605065}) between Y1 and Z1.  The function is used in \code{\link[SimMultiCorrData]{findintercorr}} and
#'     \code{\link[SimMultiCorrData]{rcorrvar}}.  This function would not ordinarily be called by the user.
#'
#' @param method the method used to generate the \code{k_cont} continuous variables.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param constants a matrix with \code{k_cont} rows, each a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or
#'     c0, c1, c2, c3, c4, c5 (if \code{method} = "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param rho_cont_pois a \code{k_cont x k_pois} matrix of target correlations among continuous and Poisson variables
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{dpois}})
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @export
#' @keywords intermediate, correlation, continuous, Fleishman, Headrick, Poisson, method 1
#' @seealso \code{\link[SimMultiCorrData]{chat_pois}}, \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{find_constants}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return a \code{k_cont x k_pois} matrix whose rows represent the \code{k_cont} continuous variables and columns represent the
#'     \code{k_pois} Poisson variables
#' @references
#' Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15): 3129-39. \doi{10.1080/00949655.2014.953534}.
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2): 104-109. \doi{10.1198/tast.2011.10090}.
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532. \doi{10.1007/BF02293811}.
#'
#' Frechet M.  Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA.  1951;14:53-77.
#'
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis, 40(4):685-711. \doi{10.1016/S0167-9473(02)00072-5}.
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3(1), 65-71. \doi{10.22237/jmasm/1083370080}.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249. \doi{10.1080/10629360600605065}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35. \doi{10.1007/BF02294317}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1): 91-102. \doi{10.1002/asmb.901}.
#'
findintercorr_cont_pois <- function(method, constants, rho_cont_pois, lam,
                                    nrand = 100000, seed = 1234) {
  Sigma_cont_pois <- matrix(1, nrow = nrow(rho_cont_pois),
                            ncol = ncol(rho_cont_pois))
  for (i in 1:nrow(rho_cont_pois)) {
    for (j in 1:ncol(rho_cont_pois)) {
      Sigma_cont_pois[i, j] <-
        rho_cont_pois[i, j]/(chat_pois(lam[j], n_unif = nrand, seed = seed) *
                               power_norm_corr(constants[i, ], method))
    }
  }
  return(Sigma_cont_pois)
}
