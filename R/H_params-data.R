#' Parameters for Examples of Constants Calculated by Headrick's Fifth-Order Polynomial Transformation
#'
#' These are the parameters for \code{\link[SimMultiCorrData]{Headrick.dist}}, which contains selected symmetrical and
#' asymmetrical theoretical densities with their associated values
#' of skewness (gamma1), standardized kurtosis (gamma2), and standardized fifth (gamma3) and
#' sixth (gamma4) cumulants.  Constants were calculated by Headrick using his fifth-order
#' polynomial transformation and given in his Table 1 (2002, p. 691-2).
#'
#' @docType data
#' @usage data(H_params)
#'
#' @format An object of class \code{"data.frame"}; Colnames are distribution names as inputs for
#' \code{\link[SimMultiCorrData]{calc_theory}}; rownames are param1, param2.
#'
#' @keywords datasets
#'
#' @references Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#' Non-normal Distributions. Computational Statistics & Data Analysis 40(4):685-711
#' (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
"H_params"
NULL
