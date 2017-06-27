#' @title Plot Theoretical Power Method Probability Density Function and Target PDF of External Data for Continuous Variables
#'
#' @description This plots the theoretical power method probability density function: \deqn{f_p(Z)(p(z)) = f_p(Z)(p(z), f_Z(z)/p'(z)),} as given in
#'     Headrick & Kowalchuk (2007, \doi{10.1080/10629360600605065}), and target
#'     pdf.  It is a parametric plot with \eqn{sigma * y + mu}, where \eqn{y = p(z)}, on the x-axis and
#'     \eqn{f_Z(z)/p'(z)} on the y-axis, where \eqn{z} is vector of \eqn{n} random standard normal numbers (generated with a seed set
#'     by user; length equal to length of external data vector).  \code{sigma} is the standard deviation and \code{mu} is the mean of the external data set.
#'     Given a vector of polynomial transformation constants, the function generates \eqn{sigma * y + mu} and calculates the theoretical
#'     probabilities using \eqn{f_p(Z)(p(z), f_Z(z)/p'(z))}.  The target distribution is also plotted given a vector
#'     of external data.  This external data is required.  The \eqn{y} values are centered and scaled to have the same mean and standard
#'     deviation as the external data.  If the user wants to only plot the power method pdf,
#'     \code{\link[SimMultiCorrData]{plot_pdf_theory}} should be used instead with \code{overlay = FALSE}.  It returns a \code{\link[ggplot2]{ggplot2}} object so the user can modify as necessary.  The graph parameters (i.e. \code{title},
#'     \code{power_color}, \code{target_color}, \code{nbins}) are \code{\link[ggplot2]{ggplot2}} parameters.  It works for valid or invalid power method pdfs.
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param method the method used to generate the continuous variable y = p(z).  "Fleishman" uses Fleishman's third-order polynomial
#'     transformation and "Polynomial" uses Headrick's fifth-order transformation.
#' @param title the title for the graph (default = "Probability Density Function")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param power_color the line color for the power method pdf (default = "dark blue")
#' @param ext_y a vector of external data (required)
#' @param target_color the histogram color for the target pdf (default = "dark green")
#' @param target_lty the line type for the target pdf (default = 2, dashed line)
#' @param seed the seed value for random number generation (default = 1234)
#' @import stats
#' @import utils
#' @import ggplot2
#' @export
#' @keywords plot, theoretical, external, pdf, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{calc_theory}},
#'     \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_line}}, \code{\link[ggplot2]{geom_density}}
#' @return A \code{\link[ggplot2]{ggplot2}} object.
#' @references Please see the references for \code{\link[SimMultiCorrData]{plot_cdf}}.
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
#' # Find constants without the sixth cumulant correction
#' # (invalid power method pdf)
#' con1 <- find_constants(method = "Polynomial", skews = stcum[3],
#'                       skurts = stcum[4], fifths = stcum[5],
#'                       sixths = stcum[6])
#'
#' # Plot invalid power method pdf with external data
#' plot_pdf_ext(c = con1$constants, method = "Polynomial",
#'              title = "Invalid Logistic PDF", ext_y = ext_y,
#'              seed = seed)
#'
#' # Find constants with the sixth cumulant correction
#' # (valid power method pdf)
#' con2 <- find_constants(method = "Polynomial", skews = stcum[3],
#'                       skurts = stcum[4], fifths = stcum[5],
#'                       sixths = stcum[6], Six = seq(1.5, 2, 0.05))
#'
#' # Plot invalid power method pdf with external data
#' plot_pdf_ext(c = con2$constants, method = "Polynomial",
#'              title = "Valid Logistic PDF", ext_y = ext_y,
#'              seed = seed)
#' }
#'
plot_pdf_ext <- function(c = NULL, method = c("Fleishman", "Polynomial"),
                         title = "Probability Density Function",
                         ylower = NULL, yupper = NULL,
                         power_color = "dark blue", ext_y = NULL,
                         target_color = "dark green", target_lty = 2,
                         seed = 1234) {
  if (is.null(ext_y)) stop("You must provide an external data set.")
  c <- as.numeric(c)
  set.seed(seed)
  n <- length(ext_y)
  z <- sort(rnorm(n, 0, 1))
  y <- numeric(n)
  phi <- numeric(n)
  fy <- numeric(n)
  sigma <- sd(ext_y)
  mu <- mean(ext_y)
  for (i in 1:length(z)) {
    phi[i] <- dnorm(z[i])
    if (method == "Fleishman") {
      y[i] <- sigma * (c[1] + c[2] * z[i] + c[3] * z[i]^2 + c[4] * z[i]^3) + mu
      fy[i] <- phi[i]/(sigma * (c[2] + 2 * c[3] * z[i] + 3 * c[4] * z[i]^2))
    }
    if (method=="Polynomial") {
      y[i] <- sigma * (c[1] + c[2] * z[i] + c[3] * z[i]^2 + c[4] * z[i]^3 +
                         c[5] * z[i]^4 + c[6] * z[i]^5) + mu
      fy[i] <- phi[i]/(sigma * (c[2] + 2 * c[3] * z[i] + 3 * c[4] * z[i]^2 +
                         4 * c[5] * z[i]^3 + 5 * c[6] * z[i]^4))
    }
  }
  if (is.null(ylower) & is.null(yupper)) {
    ylower <- min(y)
    yupper <- max(y)
  }
  data <- data.frame(y = y, fy = fy, type = as.factor(rep("sim", length(y))))
  data2 <- data.frame(y = ext_y, fy = rep(NA, length(ext_y)),
                      type = as.factor(rep("theory", length(ext_y))))
  data2 <- data.frame(rbind(data, data2))
  plot1 <- ggplot() + theme_bw() + ggtitle(title) +
    geom_line(data = data2[data2$type == "sim", ],
              aes_(x = ~y, y = ~fy, colour = ~type, lty = ~type)) +
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
                        labels = c("Power Method", "Target")) +
    scale_linetype_manual(name = "", values = c(1, target_lty),
                          labels = c("Power Method", "Target"))
  return(plot1)
}
