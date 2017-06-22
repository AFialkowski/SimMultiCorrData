#' @title Error Loop to Correct Final Correlation of Simulated Variables
#'
#' @description This function corrects the final correlation of simulated variables to be within a precision value (\code{epsilon}) of the
#'     target correlation.  It updates the pairwise intermediate MVN correlation iteratively in a loop until either the maximum error
#'     is less than epsilon or the number of iterations exceeds the maximum number set by the user (\code{maxit}).  It uses
#'     \code{\link[SimMultiCorrData]{error_vars}} to simulate the pair of variables in each iteration.  This function would not
#'     ordinarily be called directly by the user.The function is a
#'     modification of  Barbiero & Ferrari's \code{\link[GenOrd]{ordcont}} function in \code{\link[GenOrd]{GenOrd-package}}.
#'     The \code{\link[GenOrd]{ordcont}} has been modified in the following ways:
#'
#'     1) it works for continuous, ordinal (r >= 2 categories), and count variables
#'
#'     2) the initial correlation check has been removed because this intermediate correlation
#'     Sigma from \code{\link[SimMultiCorrData]{rcorrvar}} or \code{\link[SimMultiCorrData]{rcorrvar2}} has already been
#'     checked for positive-definiteness and used to generate variables (however, the pairwise correlation is checked in each
#'     iteration for positive-definiteness using the method of Higham (2002) and the \code{\link[Matrix]{nearPD}} function)
#'
#'     3) the final positive-definite check has been removed
#'
#'     4) the intermediate correlation update function was changed to accomodate more situations
#'
#'     5) a final "fail-safe" check was added at the end of the iteration loop where if the absolute
#'     error between the final and target pairwise correlation is still > 0.1, the intermediate correlation is set
#'     equal to the target correlation (if \code{extra_correct} = "TRUE"), and
#'
#'     6) allowing specifications for the sample size and the seed for reproducability.
#'
#' @param k_cat the number of ordinal (r >= 2 categories) variables
#' @param k_cont the number of continuous variables
#' @param k_pois the number of Poisson variables
#' @param k_nb the number of Negative Binomial variables
#' @param Y_cat the ordinal variables generated from \code{\link[SimMultiCorrData]{rcorrvar}} or \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @param Y the continuous (mean 0, variance 1) variables
#' @param Yb the continuous variables with desired mean and variance
#' @param Y_pois the Poisson variables
#' @param Y_nb the Negative Binomial variables
#' @param marginal a list of length equal \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param support a list of length equal \code{k_cat}; the i-th element is a vector of containing the r
#'     ordered support values; if not provided, the default is for the i-th element to be the vector 1, ..., r
#' @param method the method used to generate the continuous variables.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param means a vector of means for the continuous variables
#' @param vars a vector of variances
#' @param constants a matrix with \code{k_cont} rows, each a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or
#'     c0, c1, c2, c3, c4, c5 (if \code{method} = "Polynomial"), like that returned by
#'     \code{\link[SimMultiCorrData]{find_constants}}
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{dpois}})
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture)
#' @param n the sample size
#' @param seed the seed value for random number generation
#' @param epsilon the maximum acceptable error between the final and target correlation matrices;
#'     smaller epsilons take more time
#' @param maxit the maximum number of iterations to use to find the intermediate correlation; the
#'     correction loop stops when either the iteration number passes \code{maxit} or \code{epsilon} is reached
#' @param rho0 the target correlation matrix
#' @param Sigma the intermediate correlation matrix previously used in \code{\link[SimMultiCorrData]{rcorrvar}}
#'     or \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @param rho_calc the final correlation matrix calculated in \code{\link[SimMultiCorrData]{rcorrvar}}
#'     or \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @param extra_correct if "TRUE", a final "fail-safe" check is used at the end of the iteration loop where if the absolute
#'     error between the final and target pairwise correlation is still > 0.1, the intermediate correlation is set
#'     equal to the target correlation
#' @importFrom Matrix nearPD
#' @export
#' @keywords error, correlation
#' @seealso \code{\link[GenOrd]{ordcont}}, \code{\link[SimMultiCorrData]{rcorrvar}}, \code{\link[SimMultiCorrData]{rcorrvar2}},
#'     \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{findintercorr2}}
#'
#' @return A list with the following components:
#' @return \code{Sigma} the intermediate MVN correlation matrix resulting from the error loop
#' @return \code{rho_calc} the calculated final correlation matrix generated from Sigma
#' @return \code{Y_cat} the ordinal variables
#' @return \code{Y} the continuous (mean 0, variance 1) variables
#' @return \code{Yb} the continuous variables with desired mean and variance
#' @return \code{Y_pois} the Poisson variables
#' @return \code{Y_nb} the Negative Binomial variables
#' @return \code{niter} a matrix containing the number of iterations required for each variable pair
#' @references Ferrari PA, Barbiero A (2012). Simulating ordinal data, Multivariate Behavioral Research, 47(4): 566-589.
#'
#' Barbiero A, Ferrari PA (2015). GenOrd: Simulation of Discrete Random Variables with Given
#' Correlation Matrix and Marginal Distributions. R package version 1.4.0. \cr
#' \url{https://CRAN.R-project.org/package=GenOrd}
#'
#' Higham N (2002). Computing the nearest correlation matrix - a problem from finance; IMA Journal of Numerical Analysis 22: 329-343.
#'
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#' Non-normal Distributions. Computational Statistics & Data Analysis 40(4):685-711
#' (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35.
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3, 65-71.
error_loop <- function(k_cat, k_cont, k_pois, k_nb, Y_cat, Y, Yb, Y_pois, Y_nb,
                       marginal, support, method, means, vars, constants,
                       lam, size, prob, mu, n, seed, epsilon, maxit, rho0,
                       Sigma, rho_calc, extra_correct) {
  k <- k_cat + k_cont + k_pois + k_nb
  niter <- matrix(0, k, k)
  Sigmaold <- Sigma
  rho_calcold <- rho_calc
  for (q in 1:(k - 1)) {
    for (r in (q + 1):k) {
      if (rho0[q, r] == 0) {
        Sigma[q, r] <- 0
      } else {
        it <- 0
        while (abs(rho_calc[q, r] - rho0[q, r]) > epsilon & (it < maxit)) {
          if (rho0[q, r] * (rho0[q, r]/rho_calcold[q, r]) <= -1) {
            Sigma[q, r] <- Sigmaold[q, r] *
              (1 + 0.1 * (1 - Sigmaold[q, r]) *
                 -sign(rho0[q, r] - rho_calc[q, r]))
          }
          if (rho0[q, r] * (rho0[q, r]/rho_calcold[q, r]) >= 1) {
            Sigma[q, r] <- Sigmaold[q, r] *
              (1 + 0.1 * (1 - Sigmaold[q, r]) *
                 sign(rho0[q, r] - rho_calc[q, r]))
          }	else {
            Sigma[q, r] <- Sigmaold[q, r] * (rho0[q, r]/rho_calc[q, r])
          }
          Sigma[r, q] <- Sigma[q, r]
          Sigma2 <- matrix(c(1, Sigma[q, r], Sigma[r, q], 1),
                           nrow = 2, ncol = 2, byrow = T)
          if (min(eigen(Sigma2, symmetric = TRUE)$values) < 0) {
            Sigma2 <- nearPD(Sigma2, corr = T, keepDiag = T)$mat
            Sigma[q, r] <- Sigma2[1, 2]
            Sigma[r, q] <- Sigma[q, r]
          }
          EV <- error_vars(marginal = marginal, support = support,
                           method = method, means = means, vars = vars,
                           constants = constants, lam = lam,
                           size = size, prob = prob, mu = mu,
                           Sigma = Sigma, rho_calc = rho_calc, q = q,
                           r = r, k_cat = k_cat, k_cont = k_cont,
                           k_pois = k_pois, k_nb = k_nb, Y_cat = Y_cat,
                           Y = Y, Yb = Yb, Y_pois = Y_pois, Y_nb = Y_nb,
                           n = n, seed = seed)
          Sigma <- EV$Sigma
          rho_calc <- EV$rho_calc
          Y_cat <- EV$Y_cat
          Y <- EV$Y
          Yb <- EV$Yb
          Y_pois <- EV$Y_pois
          Y_nb <- EV$Y_nb
          Sigmaold[q, r] <- Sigma[q, r]
          Sigmaold[r, q] <- Sigmaold[q, r]
          it <- it + 1
        }
        if (extra_correct == TRUE) {
          if (it > maxit & (abs(rho_calc[q, r] - rho0[q, r])) > 0.1) {
            Sigma[q, r] <- Sigma[r, q] <- rho0[q, r]
            EV <- error_vars(marginal = marginal, support = support,
                             method = method, means = means, vars = vars,
                             constants = constants, lam = lam,
                             size = size, prob = prob, mu = mu,
                             Sigma = Sigma, rho_calc = rho_calc, q = q,
                             r = r, k_cat = k_cat, k_cont = k_cont,
                             k_pois = k_pois, k_nb = k_nb, Y_cat = Y_cat,
                             Y = Y, Yb = Yb, Y_pois = Y_pois, Y_nb = Y_nb,
                             n = n, seed = seed)
            Sigma <- EV$Sigma
            rho_calc <- EV$rho_calc
            Y_cat <- EV$Y_cat
            Y <- EV$Y
            Yb <- EV$Yb
            Y_pois <- EV$Y_pois
            Y_nb <- EV$Y_nb
          }
        }
        niter[q, r] <- it
        niter[r, q] <- it
      }
    }
  }
  return(list(Sigma = Sigma, rho_calc = rho_calc, Y_cat = Y_cat, Y = Y,
              Yb = Yb, Y_pois = Y_pois, Y_nb = Y_nb, niter = niter))
}
