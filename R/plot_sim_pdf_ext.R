#' @title Plot Simulated Probability Density Function and Target PDF of External Data
#'
#' @description This plots the pdf of continuous or count simulated data and overlays the target pdf computed from the
#'      given external data vector.  It returns a \code{\link[ggplot2]{ggplot2}} object so the user can modify as necessary.  The graph parameters
#'     (i.e. \code{title}, \code{power_color}, \code{target_color}, \code{target_lty}) are \code{\link[ggplot2]{ggplot2}} parameters.
#'     It works for valid or invalid power method pdfs.
#' @param sim_y a vector of simulated data
#' @param title the title for the graph (default = "Simulated Probability Density Function")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param power_color the histogram color for the simulated variable
#' @param ext_y a vector of external data
#' @param target_color the histogram color for the target pdf
#' @param target_lty the line type for the target pdf (default = 2, dashed line)
#' @import ggplot2
#' @export
#' @keywords plot, simulated, external, pdf, Fleishman, Headrick
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_density}}
#' @return A \code{\link[ggplot2]{ggplot2}} object.
#' @references Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
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
#'
#' Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
#'
#' @examples \dontrun{
#' # Logistic Distribution
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
#'                       Six = NULL, n = 10000, seed = seed)
#'
#' # Plot pdfs of simulated variable (invalid) and external data
#' plot_sim_pdf_ext(sim_y = Logvar1$continuous_variable,
#'              title = "Invalid Logistic Simulated PDF",
#'              ext_y = ext_y)
#'
#' # Simulate with the sixth cumulant correction
#' # (valid power method pdf)
#' Logvar2 <- nonnormvar1(method = "Polynomial", means = 0, vars = 1,
#'                       skews = stcum[3], skurts = stcum[4],
#'                       fifths = stcum[5], sixths = stcum[6],
#'                       Six = seq(1.5, 2, 0.05), n = 10000, seed = 1234)
#'
#' # Plot pdfs of simulated variable (valid) and external data
#' plot_sim_pdf_ext(sim_y = Logvar2$continuous_variable,
#'              title = "Valid Logistic Simulated PDF",
#'              ext_y = ext_y)
#' }
#'
plot_sim_pdf_ext <- function(sim_y,
                             title = "Simulated Probability Density Function",
                             ylower = NULL, yupper = NULL,
                             power_color = "dark blue", ext_y = NULL,
                             target_color = "dark green", target_lty = 2) {
  if (is.null(ext_y)) stop("You must provide an external data set.")
  if (is.null(ylower) & is.null(yupper)) {
    ylower <- min(sim_y)
    yupper <- max(sim_y)
  }
  data <- data.frame(y = sim_y, type = as.factor(rep("sim", length(sim_y))))
  data2 <- data.frame(y = ext_y,
                      type = as.factor(rep("theory", length(ext_y))))
  data2 <- data.frame(rbind(data, data2))
  plot1 <- ggplot() + theme_bw() + ggtitle(title) +
    geom_density(data = data2[data2$type == "sim", ],
                 aes_(x = ~y, colour = ~type, lty = ~type)) +
    geom_density(data = data2[data2$type == "theory", ],
                 aes_(x = ~y, colour = ~type, lty = ~type)) +
    scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
    scale_y_continuous(name = "Probability") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.position = c(0.975, 0.9), legend.justification = c(1, 1)) +
    scale_colour_manual(name = "", values = c(power_color, target_color),
                        labels = c("Simulated Variable", "Target")) +
    scale_linetype_manual(name = "", values = c(1, target_lty),
                          labels = c("Simulated Variable", "Target"))
  return(plot1)
}
