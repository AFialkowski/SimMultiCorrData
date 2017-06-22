#' Examples of Constants Calculated by Headrick's Fifth-Order Polynomial Transformation
#'
#' Selected symmetrical and asymmetrical theoretical densities with their associated values
#' of skewness (gamma1), standardized kurtosis (gamma2), and standardized fifth (gamma3) and
#' sixth (gamma4) cumulants.  Constants were calculated by Headrick using his fifth-order
#' polynomial transformation and given in his Table 1 (p. 691-2).
#'
#' @docType data
#'
#' @usage data(Headrick.dist)
#'
#' @format An object of class \code{"data.frame"}; Colnames are distribution names; rownames are
#' standardized cumulant names followed by c0, ..., c5.
#'
#' @keywords datasets
#'
#' @references Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#' Non-normal Distributions. Computational Statistics & Data Analysis 40(4):685-711
#' (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' @examples
#' z <- rnorm(10000)
#' g <- Headrick.dist$Gamma_a10b10[-c(1:4)]
#' gamma_a10b10 <- g[1] + g[2] * z + g[3] * z^2 + g[4] * z^3 + g[5] * z^4 +
#'                 g[6] * z^5
#' summary(gamma_a10b10)
"Headrick.dist"
