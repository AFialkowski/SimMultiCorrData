#' @title Plot Theoretical Power Method Cumulative Distribution Function for Continuous Variables
#'
#' @description This plots the theoretical power method cumulative distribution function:
#'     \deqn{F_p(Z)(p(z)) = F_p(Z)(p(z), F_Z(z)),} as given in Headrick & Kowalchuk (2007, \doi{10.1080/10629360600605065}).
#'     It is a parametric plot with \eqn{sigma * y + mu}, where \eqn{y = p(z)}, on the x-axis and \eqn{F_Z(z)} on the y-axis,
#'     where \eqn{z} is vector of \eqn{n} random standard normal numbers (generated with a seed set by user).  Given a vector of polynomial
#'     transformation constants, the function generates \eqn{sigma * y + mu} and calculates the theoretical cumulative probabilities
#'     using \eqn{F_p(Z)(p(z), F_Z(z))}.  If \code{calc_cprob} = TRUE, the cumulative probability up to \eqn{delta = sigma * y + mu} is
#'     calculated (see \code{\link[SimMultiCorrData]{cdf_prob}}) and the region on the plot is filled with a dashed horizontal
#'     line drawn at \eqn{F_p(Z)(delta)}.  The cumulative probability is stated on top of the line.  It returns a \code{\link[ggplot2]{ggplot2-package}} object so
#'     the user can modify as necessary.  The graph parameters (i.e. \code{title}, \code{color}, \code{fill}, \code{hline}) are
#'     \code{\link[ggplot2]{ggplot2-package}} parameters.  It works for valid or invalid power method pdfs.
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param method the method used to generate the continuous variable \eqn{y = p(z)}.  "Fleishman" uses Fleishman's third-order polynomial
#'     transformation and "Polynomial" uses Headrick's fifth-order transformation.
#' @param mu mean for the continuous variable (default = 0)
#' @param sigma standard deviation for the continuous variable (default = 1)
#' @param title the title for the graph (default = "Cumulative Distribution Function")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param calc_cprob if TRUE (default = FALSE), \code{\link[SimMultiCorrData]{cdf_prob}} is used to find the cumulative probability
#'     up to \eqn{delta = sigma * y + mu} and the region on the plot is filled with a dashed horizontal line drawn at \eqn{F_p(Z)(delta)}
#' @param delta the value \eqn{sigma * y + mu}, where \eqn{y = p(z)}, at which to evaluate the cumulative probability
#' @param color the line color for the cdf (default = "dark blue")
#' @param fill the fill color if \code{calc_cprob} = TRUE (default = "blue)
#' @param hline the dashed horizontal line color drawn at delta if \code{calc_cprob} = TRUE (default = "dark green")
#' @param n the number of random standard normal numbers to use in generating \eqn{y = p(z)} (default = 10000)
#' @param seed the seed value for random number generation (default = 1234)
#' @param text.size the size of the text displaying the cumulative probability up to \code{delta} if \code{calc_cprob} = TRUE
#' @param title.text.size the size of the plot title
#' @param axis.text.size the size of the axes text (tick labels)
#' @param axis.title.size the size of the axes titles
#' @param lower lower bound for \code{cdf_prob}
#' @param upper upper bound for \code{cdf_prob}
#' @import stats
#' @import utils
#' @import ggplot2
#' @import grid
#' @export
#' @keywords plot, theoretical, cdf, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{cdf_prob}},
#'     \code{\link[ggplot2]{ggplot2-package}}, \code{\link[ggplot2]{geom_path}}, \code{\link[ggplot2]{geom_abline}},
#'     \code{\link[ggplot2]{geom_ribbon}}
#' @return A \code{\link[ggplot2]{ggplot2-package}} object.
#' @references
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532. \doi{10.1007/BF02293811}.
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
#' Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
#'
#' @examples \dontrun{
#' # Logistic Distribution: mean = 0, sigma = 1
#'
#' # Find standardized cumulants
#' stcum <- calc_theory(Dist = "Logistic", params = c(0, 1))
#'
#' # Find constants without the sixth cumulant correction
#' # (invalid power method pdf)
#' con1 <- find_constants(method = "Polynomial", skews = stcum[3],
#'                       skurts = stcum[4], fifths = stcum[5],
#'                       sixths = stcum[6], n = 25, seed = 1234)
#'
#' # Plot cdf with cumulative probability calculated up to delta = 5
#' plot_cdf(c = con1$constants, method = "Polynomial",
#'          title = "Invalid Logistic CDF", calc_cprob = TRUE, delta = 5)
#'
#' # Find constants with the sixth cumulant correction
#' # (valid power method pdf)
#' con2 <- find_constants(method = "Polynomial", skews = stcum[3],
#'                       skurts = stcum[4], fifths = stcum[5],
#'                       sixths = stcum[6], Six = seq(1.5, 2, 0.05))
#'
#' # Plot cdf with cumulative probability calculated up to delta = 5
#' plot_cdf(c = con2$constants, method = "Polynomial",
#'          title = "Valid Logistic CDF", calc_cprob = TRUE, delta = 5)
#' }
#'
plot_cdf <- function(c = NULL, method = c("Fleishman", "Polynomial"), mu = 0,
                     sigma = 1, title = "Cumulative Distribution Function",
                     ylower = NULL, yupper = NULL, calc_cprob = FALSE,
                     delta = 5, color = "dark blue", fill = "blue",
                     hline = "dark green", n = 10000, seed = 1234,
                     text.size = 11, title.text.size = 15,
                     axis.text.size = 10, axis.title.size = 13,
                     lower = -1e06, upper = 1e06) {
  set.seed(seed)
  c <- as.numeric(c)
  z <- rnorm(n, 0, 1)
  y <- numeric(n)
  phi <- numeric(n)
  for (i in 1:length(z)) {
    if (method == "Fleishman") {
      y[i] <- sigma * (c[1] + c[2] * z[i] + c[3] * z[i]^2 + c[4] * z[i]^3) + mu
    }
    if (method == "Polynomial") {
      y[i] <- sigma * (c[1] + c[2] * z[i] + c[3] * z[i]^2 + c[4] * z[i]^3 +
                         c[5] * z[i]^4 + c[6] * z[i]^5) + mu
    }
    phi[i] <- pnorm(z[i])
  }
  if (is.null(ylower) & is.null(yupper)) {
    ylower <- min(y)
    yupper <- max(y)
  }
  data <- data.frame(y = y, phi = phi)
  data <- data[with(data, order(y)), ]
  plot1 <- ggplot(data, aes_(x = ~y, y = ~phi)) + theme_bw() + ggtitle(title) +
    geom_line(colour = color) +
    geom_hline(yintercept = 0, lty = 2, colour = "#333333") +
    geom_hline(yintercept = 1, lty = 2, colour = "#333333") +
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
    cprob <- cdf_prob(c = c, method = method, delta = delta, mu = mu,
                      sigma = sigma, lower = lower, upper = upper)
    data2 <- data.frame(y = data[data[, 1] <= delta, 1],
                        phi = data[data[, 1] <= delta, 2])
    text_one <- textGrob(paste("Cumulative probability = ",
                               round(cprob$cumulative_prob, 4), ", y = ",
                               round(delta, 4), sep = ""),
                         gp = gpar(fontsize = text.size, fontface = "bold",
                                   col = hline))
    plot1 <- plot1 +
      geom_area(data = data2, aes_(x = ~y, y = ~phi), fill = fill) +
      geom_hline(yintercept = cprob$cumulative_prob, lty = 2, colour = hline) +
      annotation_custom(text_one, xmin = 0.5 * (ylower + yupper),
                        xmax = 0.5 * (ylower + yupper), ymin = 1.03,
                        ymax = 1.03)
    return(plot1)
  }
}
