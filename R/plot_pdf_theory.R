#' @title Plot Theoretical Power Method Probability Density Function and Target PDF by Distribution Name or Function
#'
#' @description This plots the theoretical probability density function: \eqn{f_p(Z)(p(z)) = f_p(Z)(p(z), f_Z(z)/p'(z))} and target
#'     pdf (if overlay = TRUE).  It is a parametric plot with \eqn{sigma * y + mu}, where \eqn{y = p(z)}, on the x-axis and
#'     \eqn{f_Z(z)/p'(z)} on the y-axis, where \eqn{z} is vector of \eqn{n} random standard normal numbers (generated with a seed set by
#'     user).  Given a vector of polynomial
#'     transformation constants, the function generates \eqn{sigma * y + mu} and calculates the theoretical probabilities
#'     using \eqn{f_p(Z)(p(z), f_Z(z)/p'(z))}.  If \code{overlay} = TRUE, the target distribution is also plotted given either a
#'     distribution name (plus up to 3 parameters) or a pdf function \eqn{fx}.  If a target distribution is specified, \eqn{y} is
#'     scaled and then transformed so that it has the same mean and variance as the target distribution.
#'     It returns a \code{\link[ggplot2]{ggplot2}} object so the user can modify as necessary.  The graph parameters
#'     (i.e. \code{title}, \code{power_color}, \code{target_color}, \code{target_lty}) are \code{\link[ggplot2]{ggplot2}} parameters.
#'     It works for valid or invalid power method pdfs.
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param method the method used to generate the continuous variable \eqn{y = p(z)}.  "Fleishman" uses a third-order polynomial
#'     transformation and "Polynomial" uses Headrick's fifth-order transformation.
#' @param mu the desired mean for the continuous variable (used if \code{overlay = FALSE}, otherwise variable centered to have the
#'     same mean as the target distribution)
#' @param sigma the desired standard deviation for the continuous variable (used if \code{overlay = FALSE}, otherwise variable scaled
#'     to have the same standard deviation as the target distribution)
#' @param title the title for the graph (default = "Probability Density Function")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param power_color the line color for the power method pdf (default = "dark blue)
#' @param overlay if TRUE (default), the target distribution is also plotted given either a distribution name (and parameters)
#'     or pdf function fx (with bounds = ylower, yupper)
#' @param target_color the line color for the target pdf (default = "dark green")
#' @param target_lty the line type for the target pdf (default = 2, dashed line)
#' @param Dist name of the distribution. The possible values are: "Beta", "Chisq", "Exponential", "F", "Gamma", "Gaussian",
#'     "Laplace", "Logistic", "Lognormal", "Pareto", "Rayleigh", "t", "Triangular", "Uniform", "Weibull".
#'     Please refer to the documentation for each package (i.e. \code{\link[stats]{dgamma}})
#'     for information on appropriate parameter inputs.  The pareto (see \code{\link[VGAM]{dpareto}}), generalized
#'     rayleigh (see \code{\link[VGAM]{dgenray}}), and laplace (see \code{\link[VGAM]{dlaplace}}) distributions
#'     come from the \code{\link[VGAM]{VGAM}} package.  The triangular (see \code{\link[triangle]{dtriangle}}) distribution
#'     comes from the \code{\link[triangle]{triangle}} package.
#' @param params a vector of parameters (up to 3) for the desired distribution (keep NULL if \code{fx} supplied instead)
#' @param fx a pdf input as a function of x only, i.e. fx <- function(x) 0.5*(x-1)^2; must return a scalar
#'     (keep NULL if \code{Dist} supplied instead)
#' @param lower the lower support bound for \code{fx}
#' @param upper the upper support bound for \code{fx}
#' @param n the number of random standard normal numbers to use in generating \eqn{y = p(z)} (default = 100)
#' @param seed the seed value for random number generation (default = 1234)
#' @import ggplot2
#' @import stats
#' @import utils
#' @importFrom VGAM dpareto rpareto dgenray rgenray dlaplace rlaplace
#' @importFrom triangle dtriangle rtriangle
#' @export
#' @keywords plot, theoretical, pdf, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{calc_theory}},
#'     \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_line}}
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
#' # Find standardized cumulants
#' stcum <- calc_theory(Dist = "Logistic", params = c(0, 1))
#'
#' # Find constants without the sixth cumulant correction
#' # (invalid power method pdf)
#' con1 <- find_constants(method = "Polynomial", skews = stcum[3],
#'                       skurts = stcum[4], fifths = stcum[5],
#'                       sixths = stcum[6])
#'
#' # Plot invalid power method pdf with theoretical pdf overlayed
#' plot_pdf_theory(c = con1$constants, method = "Polynomial",
#'          title = "Invalid Logistic PDF", overlay = TRUE,
#'          Dist = "Logistic", params = c(0, 1))
#'
#' # Find constants with the sixth cumulant correction
#' # (valid power method pdf)
#' con2 <- find_constants(method = "Polynomial", skews = stcum[3],
#'                       skurts = stcum[4], fifths = stcum[5],
#'                       sixths = stcum[6], Six = seq(1.5, 2, 0.05))
#'
#' # Plot valid power method pdf with theoretical pdf overlayed
#' plot_pdf_theory(c = con2$constants, method = "Polynomial",
#'          title = "Valid Logistic PDF", overlay = TRUE,
#'          Dist = "Logistic", params = c(0, 1))
#' }
#'
plot_pdf_theory <- function(c = NULL, method = c("Fleishman", "Polynomial"),
                            mu = 0, sigma = 1,
                            title = "Probability Density Function",
                            ylower = NULL, yupper = NULL,
                            power_color = "dark blue", overlay = TRUE,
                            target_color = "dark green", target_lty = 2,
                            Dist = c("Beta", "Chisq", "Exponential", "F",
                                     "Gamma", "Gaussian", "Laplace",
                                     "Logistic", "Lognormal", "Pareto",
                                     "Rayleigh", "t", "Triangular", "Uniform",
                                     "Weibull"), params = NULL,
                            fx = NULL, lower = NULL, upper = NULL,
                            n = 100, seed = 1234) {
  if (overlay == FALSE) {
    set.seed(seed)
    c <- as.numeric(c)
    z <- rnorm(n, 0, 1)
    z <- sort(z)
    y <- numeric(n)
    phi <- numeric(n)
    fy <- numeric(n)
    for (i in 1:length(z)) {
      phi[i] <- dnorm(z[i])
      if (method == "Fleishman") {
        y[i] <- sigma * (c[1] + c[2] * z[i] + c[3] * z[i]^2 + c[4] * z[i]^3) + mu
        fy[i] <- phi[i]/(sigma * (c[2] + 2 * c[3] * z[i] + 3 * c[4] * z[i]^2))
      }
      if (method == "Polynomial") {
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
    plot1 <- ggplot(data, aes_(x = ~y, y = ~fy, colour = ~type)) + theme_bw() +
      ggtitle(title) + geom_line() +
      scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
      scale_y_continuous(name = "Probability") +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 13),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 13),
            legend.text = element_text(size = 10),
            legend.position = c(0.975, 0.9), legend.justification=c(1, 1)) +
      scale_colour_manual(name = "", values = c(power_color),
                          labels = "Power Method")
    return(plot1)
  }
  if (overlay == TRUE) {
    if (!is.null(fx)) {
      theory <- calc_theory(fx = fx, lower = lower, upper = upper)
    }
    if (is.null(fx)) {
      theory <- calc_theory(Dist = Dist, params = params)
    }
    mu <- theory[1]
    sigma <- theory[2]
    set.seed(seed)
    c <- as.numeric(c)
    z <- rnorm(n, 0, 1)
    z <- sort(z)
    y <- numeric(n)
    phi <- numeric(n)
    fy <- numeric(n)
    for (i in 1:length(z)) {
      phi[i] <- dnorm(z[i])
      if (method == "Fleishman") {
        y[i] <- sigma * (c[1] + c[2] * z[i] + c[3] * z[i]^2
                         + c[4] * z[i]^3) + mu
        fy[i] <- phi[i]/(sigma * (c[2] + 2 * c[3] * z[i] + 3 * c[4] * z[i]^2))
      }
      if (method == "Polynomial") {
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
    x <- y
    y_fx <- numeric(n)
    if (!is.null(fx)) {
      for (j in 1:n) {
        y_fx[j] <- fx(x[j])
      }
    }
    if (is.null(fx)) {
      D <- data.frame(Dist = c("Beta", "Chisq", "Exponential", "F", "Gamma",
                               "Gaussian", "Laplace", "Logistic", "Lognormal",
                               "Pareto", "Rayleigh", "t", "Triangular",
                               "Uniform", "Weibull"),
                      pdf = c("dbeta", "dchisq", "dexp", "df", "dgamma",
                              "dnorm", "dlaplace", "dlogis", "dlnorm",
                              "dpareto", "dgenray", "dt", "dtriangle",
                              "dunif", "dweibull"),
                      fx = c("rbeta", "rchisq", "rexp", "rf", "rgamma",
                             "rnorm", "rlaplace", "rlogis", "rlnorm",
                             "rpareto", "rgenray", "rt", "rtriangle",
                             "runif", "rweibull"),
                      Lower = as.numeric(c(0, 0, 0, 0, 0, -Inf, -Inf, -Inf,
                                           0, params[1], 0, -Inf, params[1],
                                           params[1], 0)),
                      Upper = as.numeric(c(1, Inf, Inf, Inf, Inf, Inf, Inf,
                                           Inf, Inf, Inf, Inf, Inf, params[2],
                                           params[2], Inf)))
      i <- match(Dist, D$Dist)
      p <- as.character(D$pdf[i])
      for (j in 1:n) {
        if (length(params) == 1) y_fx[j] <- get(p)(x[j], params[1])
        if (length(params) == 2) y_fx[j] <- get(p)(x[j], params[1], params[2])
        if (length(params) == 3) y_fx[j] <- get(p)(x[j], params[1], params[2],
                                                   params[3])
      }
    }
    data <- data.frame(y = y, fy = fy, type = as.factor(rep("sim", length(y))))
    data2 <- data.frame(y = x, fy = y_fx, type = as.factor(rep("theory",
                                                               length(y_fx))))
    data2 <- data.frame(rbind(data, data2))
    plot1 <- ggplot() + theme_bw() + ggtitle(title) +
      geom_line(data = data2[data2$type == "sim", ],
                aes_(x = ~y, y = ~fy, colour = ~type, lty = ~type)) +
      geom_line(data = data2[data2$type == "theory", ],
                aes_(x = ~y, y = ~fy, colour = ~type, lty = ~type)) +
      scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
      scale_y_continuous(name = "Probability",
                         limits = c(min(data2$fy), max(data2$fy))) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 13),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 13),
            legend.text = element_text(size = 10),
            legend.position = c(0.975, 0.9), legend.justification = c(1, 1)) +
      scale_colour_manual(name = "", values = c(power_color, target_color),
                          labels = c("Power Method", "Target")) +
      scale_linetype_manual(name ="", values = c(1, target_lty),
                            labels = c("Power Method", "Target"))
    return(plot1)
  }
}
