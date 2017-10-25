#' @title Plot Simulated (Empirical) Cumulative Distribution Function for Continuous, Ordinal, or Count Variables
#'
#' @description This plots the cumulative distribution function of simulated continuous, ordinal, or count data using the empirical cdf
#'     \eqn{Fn} (see \code{\link[ggplot2]{stat_ecdf}}).
#'     \eqn{Fn} is a step function with jumps \eqn{i/n} at observation values, where \eqn{i} is the number of tied observations at that
#'     value. Missing values are
#'     ignored.  For observations \eqn{y = (y1, y2, ..., yn)}, \eqn{Fn} is the fraction of observations less or equal to \eqn{t}, i.e.,
#'     \eqn{Fn(t) = sum[yi <= t]/n}.  If \code{calc_cprob} = TRUE and the variable is \emph{continuous}, the cumulative probability up to
#'     \eqn{y = delta} is calculated (see \code{\link[SimMultiCorrData]{sim_cdf_prob}}) and the region on the plot is filled with a
#'     dashed horizontal line drawn at Fn(delta).  The cumulative probability is stated on top of the line.
#'     This fill option does not work for ordinal or count variables.  The function returns a
#'     \code{\link[ggplot2]{ggplot2}} object so the user can modify as necessary.
#'     The graph parameters (i.e. \code{title}, \code{color}, \code{fill}, \code{hline}) are \code{\link[ggplot2]{ggplot2}} parameters.
#'     It works for valid or invalid power method pdfs.
#' @param sim_y a vector of simulated data
#' @param title the title for the graph (default = "Empirical Cumulative Distribution Function")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param calc_cprob if TRUE (default = FALSE) and \code{sim_y} is continuous, \code{\link[SimMultiCorrData]{sim_cdf_prob}} is used to find the empirical cumulative probability
#'     up to y = delta and the region on the plot is filled with a dashed horizontal line drawn at \eqn{Fn(delta)}
#' @param delta the value y at which to evaluate the cumulative probability (default = 5)
#' @param color the line color for the cdf (default = "dark blue")
#' @param fill the fill color if \code{calc_cprob} = TRUE (default = "blue)
#' @param hline the dashed horizontal line color drawn at \code{delta} if \code{calc_cprob} = TRUE (default = "dark green")
#' @param text.size the size of the text displaying the cumulative probability up to \code{delta} if \code{calc_cprob} = TRUE
#' @param title.text.size the size of the plot title
#' @param axis.text.size the size of the axes text (tick labels)
#' @param axis.title.size the size of the axes titles
#' @import ggplot2
#' @import grid
#' @export
#' @keywords plot, simulated, empirical, cdf
#' @seealso \code{\link[stats]{ecdf}}, \code{\link[SimMultiCorrData]{sim_cdf_prob}}, \code{\link[ggplot2]{ggplot}},
#'     \code{\link[ggplot2]{stat_ecdf}}, \code{\link[ggplot2]{geom_hline}}, \code{\link[ggplot2]{geom_area}}
#' @return A \code{\link[ggplot2]{ggplot2}} object.
#' @references Please see the references for \code{\link[SimMultiCorrData]{plot_cdf}}.
#'
#' Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
#'
#' @examples \dontrun{
#' # Logistic Distribution: mean = 0, variance = 1
#' seed = 1234
#'
#' # Find standardized cumulants
#' stcum <- calc_theory(Dist = "Logistic", params = c(0, 1))
#'
#' # Simulate without the sixth cumulant correction
#' # (invalid power method pdf)
#' Logvar1 <- nonnormvar1(method = "Polynomial", means = 0, vars = 1,
#'                       skews = stcum[3], skurts = stcum[4],
#'                       fifths = stcum[5], sixths = stcum[6], seed = seed)
#'
#' # Plot cdf with cumulative probability calculated up to delta = 5
#' plot_sim_cdf(sim_y = Logvar1$continuous_variable,
#'              title = "Invalid Logistic Empirical CDF",
#'              calc_cprob = TRUE, delta = 5)
#'
#' # Simulate with the sixth cumulant correction
#' # (valid power method pdf)
#' Logvar2 <- nonnormvar1(method = "Polynomial", means = 0, vars = 1,
#'                       skews = stcum[3], skurts = stcum[4],
#'                       fifths = stcum[5], sixths = stcum[6],
#'                       Six = seq(1.5, 2, 0.05), seed = seed)
#'
#' # Plot cdf with cumulative probability calculated up to delta = 5
#' plot_sim_cdf(sim_y = Logvar2$continuous_variable,
#'              title = "Valid Logistic Empirical CDF",
#'              calc_cprob = TRUE, delta = 5)
#'
#' # Simulate one binary and one ordinal variable (4 categories) with
#' # correlation 0.3
#' Ordvars = rcorrvar(k_cat = 2, marginal = list(0.4, c(0.2, 0.5, 0.7)),
#'                    rho = matrix(c(1, 0.3, 0.3, 1), 2, 2), seed = seed)
#'
#' # Plot cdf of 2nd variable
#' plot_sim_cdf(Ordvars$ordinal_variables[, 2])
#'
#' }
#'
plot_sim_cdf <- function(sim_y,
                         title = "Empirical Cumulative Distribution Function",
                         ylower = NULL, yupper = NULL, calc_cprob = FALSE,
                         delta = 5, color = "dark blue", fill = "blue",
                         hline = "dark green", text.size = 11,
                         title.text.size = 15, axis.text.size = 10,
                         axis.title.size = 13) {
  if (is.null(ylower) & is.null(yupper)) {
    ylower <- min(sim_y)
    yupper <- max(sim_y)
  }
  data <- data.frame(y = sim_y)
  plot1 <- ggplot(data,  aes_(~y)) + stat_ecdf(geom = "step", colour = color) +
    theme_bw() + geom_hline(yintercept = 0, lty = 2, colour = "#333333") +
    geom_hline(yintercept = 1, lty = 2, colour = "#333333") + ggtitle(title) +
    scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
    scale_y_continuous(name = "Cumulative Probability") +
    theme(plot.title = element_text(size = title.text.size, face = "bold",
                                    hjust = 0.5),
          legend.position = "none",
          axis.text.x  = element_text(size = axis.text.size),
          axis.title.x = element_text(size = axis.title.size),
          axis.text.y  = element_text(size = axis.text.size),
          axis.title.y = element_text(size = axis.title.size))
  if (calc_cprob == FALSE) return(plot1)
  if (calc_cprob == TRUE) {
    cprob <- sim_cdf_prob(sim_y = sim_y, delta = delta)
    data2 <- data.frame(y = data[data[, 1] <= delta, ])
    data2$cum_prob <- cprob$Fn(as.numeric(data2[, 1]))
    text_one <- textGrob(paste("Cumulative probability = ",
                               round(cprob$cumulative_prob, 4), ", y = ",
                               round(delta, 4), sep = ""),
                         gp = gpar(fontsize = text.size, fontface = "bold",
                                   col = hline))
    plot1 <- plot1 +
      geom_area(data = data2, aes_(x = ~y, y = ~cum_prob), fill = fill) +
      geom_hline(yintercept = cprob$cumulative_prob, lty = 2, colour = hline) +
      annotation_custom(text_one, xmin = 0.5 * (ylower + yupper),
                        xmax = 0.5 * (ylower + yupper), ymin = 1.03,
                        ymax = 1.03)
    return(plot1)
  }
}
