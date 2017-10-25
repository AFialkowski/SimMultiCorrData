#' @title Calculate Intermediate MVN Correlation for Continuous - Poisson Variables: Correlation Method 2
#'
#' @description This function calculates a \code{k_cont x k_pois} intermediate matrix of correlations for the \code{k_cont} continuous and
#'     \code{k_pois} Poisson variables. It extends the methods of Demirtas et al. (2012, \doi{10.1002/sim.5362}) and
#'     Barbiero & Ferrari (2015, \doi{10.1002/asmb.2072}) by:
#'
#'     1) including non-normal continuous and count variables
#'
#'     2) allowing the continuous variables to be generated via Fleishman's third-order or Headrick's fifth-order transformation, and
#'
#'     3) since the count variables are treated as ordinal, using the point-polyserial and polyserial correlations to calculate the
#'     intermediate correlations (similar to \code{\link[SimMultiCorrData]{findintercorr_cont_cat}}).
#'
#'     Here, the intermediate correlation between Z1 and Z2 (where Z1
#'     is the standard normal variable transformed using Headrick's fifth-order or Fleishman's third-order method to produce a
#'     continuous variable Y1, and Z2 is the standard normal variable used to generate a Poisson variable via the inverse cdf method)
#'     is calculated by dividing the target correlation by a correction factor.  The correction factor is the product of the
#'     point-polyserial correlation between Y2 and Z2 (described in Olsson et al., 1982, \doi{10.1007/BF02294164})
#'     and the power method correlation (described in Headrick & Kowalchuk, 2007, \doi{10.1080/10629360600605065}) between Y1 and Z1.
#'     After the maximum support value has been found using
#'     \code{\link[SimMultiCorrData]{max_count_support}}, the point-polyserial correlation is given by:
#'     \deqn{\rho_{y2,z2} = (1/\sigma_{y2})\sum_{j = 1}^{r-1} \phi(\tau_{j})(y2_{j+1} - y2_{j})} where
#'     \deqn{\phi(\tau) = (2\pi)^{-1/2}*exp(-\tau^2/2)}  Here, \eqn{y_{j}} is the j-th support
#'     value and \eqn{\tau_{j}} is \eqn{\Phi^{-1}(\sum_{i=1}^{j} Pr(Y = y_{i}))}.  The power method correlation is given by:
#'     \deqn{\rho_{y1,z1} = c1 + 3c3 + 15c5}, where c5 = 0 if \code{method} = "Fleishman".  The function is used in
#'     \code{\link[SimMultiCorrData]{findintercorr2}} and \code{\link[SimMultiCorrData]{rcorrvar2}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param method the method used to generate the \code{k_cont} continuous variables.  "Fleishman" uses Fleishman's third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param constants a matrix with \code{k_cont} rows, each a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or
#'     c0, c1, c2, c3, c4, c5 (if \code{method} = "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param rho_cont_pois a \code{k_cont x k_pois} matrix of target correlations among continuous and Poisson variables
#' @param pois_marg a list of length equal to \code{k_pois}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1);
#'     this is created within \code{\link[SimMultiCorrData]{findintercorr2}} and \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @param pois_support a list of length equal to \code{k_pois}; the i-th element is a vector of containing the r
#'     ordered support values, with a minimum of 0 and maximum determined via \code{\link[SimMultiCorrData]{max_count_support}}
#' @export
#' @keywords intermediate, correlation, continuous, Fleishman, Headrick, Poisson, method 2
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{findintercorr2}}, \code{\link[SimMultiCorrData]{rcorrvar2}}
#' @return a \code{k_cont x k_pois} matrix whose rows represent the \code{k_cont} continuous variables and columns represent the
#'     \code{k_pois} Poisson variables
#' @references Please see additional references in \code{\link[SimMultiCorrData]{findintercorr_cont_cat}}.
#'
#' Barbiero A & Ferrari PA (2015). Simulation of correlated Poisson variables. Applied Stochastic Models in
#'     Business and Industry, 31: 669-80. \doi{10.1002/asmb.2072}.
#'
findintercorr_cont_pois2 <- function(method, constants, rho_cont_pois,
                                     pois_marg, pois_support) {
  Sigma_cont_pois <- matrix(1, nrow = nrow(rho_cont_pois),
                            ncol = ncol(rho_cont_pois))
  for (i in 1:nrow(rho_cont_pois)) {
    for (j in 1:ncol(rho_cont_pois)) {
      Sigma_cont_pois[i, j] <-
      (rho_cont_pois[i, j] * sqrt(var_cat(marginal = pois_marg[[j]],
      support = pois_support[[j]])))/(denom_corr_cat(marginal = pois_marg[[j]],
      support = pois_support[[j]]) * power_norm_corr(constants[i, ], method))
    }
  }
  return(Sigma_cont_pois)
}
