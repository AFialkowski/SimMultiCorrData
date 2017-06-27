#' @title Calculate Intermediate MVN Correlation to Generate Variables Treated as Ordinal
#'
#' @description This function calculates the intermediate MVN correlation needed to generate a variable described by
#'     a discrete marginal distribution and associated finite support.  This includes ordinal (r >= 2 categories) variables
#'     or variables that are treated as ordinal (i.e. count variables in the Barbiero & Ferrari, 2015 method used in
#'     \code{\link[SimMultiCorrData]{rcorrvar2}}, \doi{10.1002/asmb.2072}).  The function is a modification of  Barbiero & Ferrari's
#'     \code{\link[GenOrd]{ordcont}}
#'     function in \code{\link[GenOrd]{GenOrd-package}}.  It works by setting the intermediate MVN correlation equal to the target
#'     correlation and updating each intermediate pairwise correlation until the final pairwise correlation is within epsilon of the
#'     target correlation or the maximum number of iterations has been reached.  This function uses \code{\link[GenOrd]{contord}}
#'     to calculate the ordinal correlation obtained from discretizing the normal variables generated from the intermediate
#'     correlation matrix.  The \code{\link[GenOrd]{ordcont}} has been modified in the following ways:
#'
#'     1) the initial correlation check has been removed because it is assumed the user has done this before simulation using
#'     \code{\link[SimMultiCorrData]{valid_corr}} or \code{\link[SimMultiCorrData]{valid_corr2}}
#'
#'     2) the final positive-definite check has been removed
#'
#'     3) the intermediate correlation update function was changed to accomodate more situations, and
#'
#'     4) a final "fail-safe" check was added at the end of the iteration loop where if the absolute
#'     error between the final and target pairwise correlation is still > 0.1, the intermediate correlation is set
#'     equal to the target correlation.
#'
#'     This function would not ordinarily be called by the user.  Note that this will return a matrix that is NOT positive-definite
#'     because this is corrected for in the
#'     simulation functions \code{\link[SimMultiCorrData]{rcorrvar}} and \code{\link[SimMultiCorrData]{rcorrvar2}}
#'     using the method of Higham (2002) and the \code{\link[Matrix]{nearPD}} function.
#' @param marginal a list of length equal to the number of variables; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param rho the target correlation matrix
#' @param support a list of length equal to the number of variables; the i-th element is a vector of containing the r
#'     ordered support values; if not provided (i.e. support = list()), the default is for the i-th element to be the vector 1, ..., r
#' @param epsilon the maximum acceptable error between the final and target correlation matrices (default = 0.001);
#'     smaller epsilons take more time
#' @param maxit the maximum number of iterations to use (default = 1000) to find the intermediate correlation; the
#'     correction loop stops when either the iteration number passes maxit or epsilon is reached
#' @import GenOrd
#' @export
#' @keywords intermediate, correlation, ordinal, count
#' @seealso \code{\link[GenOrd]{ordcont}}, \code{\link[SimMultiCorrData]{rcorrvar}}, \code{\link[SimMultiCorrData]{rcorrvar2}},
#'          \code{\link[SimMultiCorrData]{findintercorr}}, \code{\link[SimMultiCorrData]{findintercorr2}}
#' @return A list with the following components:
#' @return \code{SigmaC} the intermediate MVN correlation matrix
#' @return \code{rho0} the calculated final correlation matrix generated from \code{SigmaC}
#' @return \code{rho} the target final correlation matrix
#' @return \code{niter} a matrix containing the number of iterations required for each variable pair
#' @return \code{maxerr} the maximum final error between the final and target correlation matrices
#' @references
#' Barbiero A, Ferrari PA (2015). Simulation of correlated Poisson variables. Applied Stochastic Models
#'     in Business and Industry, 31: 669-80. \doi{10.1002/asmb.2072}.
#'
#' Barbiero A, Ferrari PA (2015). GenOrd: Simulation of Discrete Random Variables with Given
#'     Correlation Matrix and Marginal Distributions. R package version 1.4.0.
#'     \url{https://CRAN.R-project.org/package=GenOrd}
#'
#' Ferrari PA, Barbiero A (2012). Simulating ordinal data, Multivariate Behavioral Research, 47(4): 566-589. \doi{10.1080/00273171.2012.692630}.
#'
ordnorm <- function(marginal, rho, support = list(), epsilon = 0.001,
                    maxit = 1000) {
  if (!all(unlist(lapply(marginal,
                         function(x) (sort(x) == x & min(x) > 0 &
                                      max(x) < 1))))) {
    stop("Error in given marginal distributions!")
  }
  if (!isSymmetric(rho) | min(eigen(rho, symmetric = TRUE)$values) < 0 |
      !all(diag(rho) == 1)) {
    stop("Correlation matrix not valid!")
  }
  len <- length(support)
  k_cat <- length(marginal)
  niter <- matrix(0, k_cat, k_cat)
  kj <- numeric(k_cat)
  for (i in 1:k_cat) {
    kj[i] <- length(marginal[[i]]) + 1
    if (len == 0) support[[i]] <- 1:kj[i]
  }
  rho0 <- rho
  rhoord <- contord(marginal, rho, support, Spearman = FALSE)
  rhoold <- rho
  rhoordold <- rhoord
  for (q in 1:(k_cat - 1)) {
    for (r in (q + 1):k_cat) {
      if (rho0[q, r] == 0) {
        rho[q, r] <- 0
      }
      else {
        it <- 0
        while ((abs(rhoord[q, r] - rho0[q, r]) > epsilon) & it < maxit) {
          if (rho0[q, r] * (rho0[q, r]/rhoordold[q, r]) <= -1) {
            rho[q, r] <- rhoold[q, r] * (1 + 0.1 * (1 - rhoold[q, r]) *
                                           -sign(rho0[q, r] - rhoord[q, r]))
          }
          if (rho0[q, r] * (rho0[q, r]/rhoordold[q, r]) >= 1) {
            rho[q, r] <- rhoold[q, r] * (1 + 0.1 * (1 - rhoold[q, r]) *
                                           sign(rho0[q, r] - rhoord[q, r]))
          }	else {
            rho[q, r] <- rhoold[q, r] * (rho0[q, r]/rhoord[q, r])
          }
          rho[r, q] <- rho[q, r]
          rhoord[r, q] <- contord(list(marginal[[q]], marginal[[r]]),
                                  matrix(c(1, rho[q, r], rho[q, r], 1), 2, 2),
                                  list(support[[q]], support[[r]]),
                                  Spearman = FALSE)[2]
          rhoord[q, r] <- rhoord[r, q]
          rhoold[q, r] <- rho[q, r]
          rhoold[r, q] <- rhoold[q, r]
          it <- it + 1
        }
        if (it > maxit & abs(rhoord[q, r] - rho0[q, r]) > 0.1) {
          rho[q, r] <- rho[r, q] <- rho0[q, r]
        }
        niter[q, r] <- it
        niter[r, q] <- it
      }
    }
  }
  emax <- max(abs(rhoord - rho0))
  list(SigmaC = rho, rhoO = rhoord, rho = rho0, niter = niter, maxerr = emax)
}
