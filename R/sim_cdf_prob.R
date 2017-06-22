#' @title Calculate Simulated (Empirical) Cumulative Probability
#'
#' @description This function calculates a cumulative probability using simulated data and
#'     Martin Maechler's \code{\link[stats]{ecdf}} function.  Fn is a step function with jumps i/n at observation
#'     values, where i is the number of tied observations at that value. Missing values are ignored. For
#'     observations y= (y1,y2, ... yn), Fn is the fraction of observations less or equal to t, i.e.,
#'     \eqn{Fn(t) = sum[yi <= t]/n}.
#' @param sim_y a vector of simulated data
#' @param delta the value y at which to evaluate the cumulative probability
#' @import stats
#' @import utils
#' @export
#' @keywords simulated, empirical, statistics, cumulative, probability
#' @seealso \code{\link[stats]{ecdf}}, \code{\link[SimMultiCorrData]{plot_sim_cdf}}
#' @return A list with components:
#' @return \code{cumulative_prob} the empirical cumulative probability up to delta
#' @return \code{Fn} the empirical distribution function
#'
#' @examples
#' # Beta(a = 4, b = 2) Distribution:
#' x <- rbeta(10000, 4, 2)
#' sim_cdf_prob(x, delta = 0.5)
#'
sim_cdf_prob <- function(sim_y, delta = 0.5) {
  Fn <- ecdf(sim_y)
  cum_prob <- Fn(delta)
  return(list(cumulative_prob = cum_prob, Fn = Fn))
}
