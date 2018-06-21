#' @title Plot Simulated Data and Target External Data for Continuous or Count Variables
#'
#' @description This plots simulated continuous or count data and overlays external data, both as histograms.
#'     The external data is a required input.  The simulated data is centered and scaled to have the same mean and variance as the external
#'     data set.  If the user wants to only plot simulated data, \code{\link[SimMultiCorrData]{plot_sim_theory}} should be used instead
#'     with \code{overlay = FALSE}.
#'     It returns a \code{\link[ggplot2]{ggplot2-package}} object so the user can modify as necessary.
#'     The graph parameters (i.e. \code{title}, \code{power_color}, \code{target_color}, \code{nbins}) are
#'     \code{\link[ggplot2]{ggplot2-package}} parameters.  It works for valid or invalid power method pdfs.
#' @param sim_y a vector of simulated data
#' @param title the title for the graph (default = "Simulated Data Values")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param power_color the histogram fill color for the simulated variable (default = "dark blue")
#' @param ext_y a vector of external data (required)
#' @param target_color the histogram fill color for the target data (default = "dark green")
#' @param nbins the number of bins to use in generating the histograms (default = 100)
#' @param legend.position the position of the legend
#' @param legend.justification the justification of the legend
#' @param legend.text.size the size of the legend labels
#' @param title.text.size the size of the plot title
#' @param axis.text.size the size of the axes text (tick labels)
#' @param axis.title.size the size of the axes titles
#' @import ggplot2
#' @export
#' @keywords plot, simulated, external, Fleishman, Headrick
#' @seealso \code{\link[ggplot2]{ggplot2-package}}, \code{\link[ggplot2]{geom_histogram}}
#' @return A \code{\link[ggplot2]{ggplot2-package}} object.
#' @references Please see the references for \code{\link[SimMultiCorrData]{plot_cdf}}.
#'
#' Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
#'
#' @examples \dontrun{
#' # Logistic Distribution: mean = 0, variance = 1
#'
#' seed = 1234
#'
#' # Simulate "external" data set
#' set.seed(seed)
#' ext_y <- rlogis(10000)
#'
#' # Find standardized cumulants
#' stcum <- calc_theory(Dist = "Logistic", params = c(0, 1))
#'
#' # Simulate without the sixth cumulant correction
#' # (invalid power method pdf)
#' Logvar1 <- nonnormvar1(method = "Polynomial", means = 0, vars = 1,
#'                       skews = stcum[3], skurts = stcum[4],
#'                       fifths = stcum[5], sixths = stcum[6],
#'                       n = 10000, seed = seed)
#'
#' # Plot simulated variable and external data
#' plot_sim_ext(sim_y = Logvar1$continuous_variable,
#'              title = "Invalid Logistic Simulated Data Values",
#'              ext_y = ext_y)
#'
#' # Simulate with the sixth cumulant correction
#' # (valid power method pdf)
#' Logvar2 <- nonnormvar1(method = "Polynomial", means = 0, vars = 1,
#'                       skews = stcum[3], skurts = stcum[4],
#'                       fifths = stcum[5], sixths = stcum[6],
#'                       Six = seq(1.5, 2, 0.05), n = 10000, seed = seed)
#'
#' # Plot simulated variable and external data
#' plot_sim_ext(sim_y = Logvar2$continuous_variable,
#'              title = "Valid Logistic Simulated Data Values",
#'              ext_y = ext_y)
#'
#' # Simulate 2 Poisson distributions (means = 10, 15) and correlation 0.3
#' # using Method 1
#' Pvars <- rcorrvar(k_pois = 2, lam = c(10, 15),
#'                   rho = matrix(c(1, 0.3, 0.3, 1), 2, 2), seed = seed)
#'
#' # Simulate "external" data set
#' set.seed(seed)
#' ext_y <- rpois(10000, 10)
#'
#' # Plot 1st simulated variable and external data
#' plot_sim_ext(sim_y = Pvars$Poisson_variable[, 1], ext_y = ext_y)
#'
#' }
#'
plot_sim_ext <- function(sim_y, title = "Simulated Data Values",
                         ylower = NULL, yupper = NULL,
                         power_color = "dark blue", ext_y = NULL,
                         target_color = "dark green", nbins = 100,
                         legend.position = c(0.975, 0.9),
                         legend.justification = c(1, 1),
                         legend.text.size = 10, title.text.size = 15,
                         axis.text.size = 10, axis.title.size = 13) {
  if (is.null(ext_y)) stop("You must provide an external data set.")
  sim_y <- sd(ext_y) * scale(sim_y) + mean(ext_y)
  if (is.null(ylower) & is.null(yupper)) {
    ylower <- min(sim_y)
    yupper <- max(sim_y)
  }
  data <- data.frame(x = 1:length(sim_y), y = sort(sim_y),
                     type = as.factor(rep("sim", length(sim_y))))
  data2 <- data.frame(x = 1:length(ext_y), y = sort(ext_y),
                      type = as.factor(rep("theory", length(ext_y))))
  data2 <- data.frame(rbind(data, data2))
  plot1 <- ggplot() + theme_bw() + ggtitle(title) +
    geom_histogram(data = data2[data2$type == "sim", ],
                   aes_(~y, fill = ~type), bins = nbins) +
    geom_histogram(data = data2[data2$type == "theory", ],
                   aes_(~y, fill = ~type), bins = nbins) +
    scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
    theme(plot.title = element_text(size = title.text.size, face = "bold",
                                    hjust = 0.5),
          axis.text.x = element_text(size = axis.text.size),
          axis.title.x = element_text(size = axis.title.size),
          axis.text.y = element_text(size = axis.text.size),
          axis.title.y = element_text(size = axis.title.size),
          legend.text = element_text(size = legend.text.size),
          legend.position = legend.position,
          legend.justification = legend.justification) +
    scale_fill_manual(name = "", values = c("sim" = power_color,
                                            "theory" = target_color),
                      labels = c("Simulated Variable", "Target Variable"))
  return(plot1)
}
