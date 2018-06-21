#' @title Plot Theoretical Power Method Probability Density Function and Target PDF by Distribution Name or Function for Continuous Variables
#'
#' @description This plots the theoretical power method probability density function: \deqn{f_p(Z)(p(z)) = f_p(Z)(p(z), f_Z(z)/p'(z)),} as given
#'     in Headrick & Kowalchuk (2007, \doi{10.1080/10629360600605065}), and target
#'     pdf (if overlay = TRUE).  It is a parametric plot with \eqn{sigma * y + mu}, where \eqn{y = p(z)}, on the x-axis and
#'     \eqn{f_Z(z)/p'(z)} on the y-axis, where \eqn{z} is vector of \eqn{n} random standard normal numbers (generated with a seed set by
#'     user).  Given a vector of polynomial
#'     transformation constants, the function generates \eqn{sigma * y + mu} and calculates the theoretical probabilities
#'     using \eqn{f_p(Z)(p(z), f_Z(z)/p'(z))}.  If \code{overlay} = TRUE, the target distribution is also plotted given either a
#'     distribution name (plus up to 4 parameters) or a pdf function \eqn{fx}.  If a target distribution is specified, \eqn{y} is
#'     scaled and then transformed so that it has the same mean and variance as the target distribution.
#'     It returns a \code{\link[ggplot2]{ggplot2-package}} object so the user can modify as necessary.  The graph parameters
#'     (i.e. \code{title}, \code{power_color}, \code{target_color}, \code{target_lty}) are \code{\link[ggplot2]{ggplot2-package}} parameters.
#'     It works for valid or invalid power method pdfs.
#' @param c a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or c0, c1, c2, c3, c4, c5 (if \code{method} =
#'     "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param method the method used to generate the continuous variable \eqn{y = p(z)}.  "Fleishman" uses Fleishman's third-order polynomial
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
#' @param Dist name of the distribution. The possible values are: "Benini", "Beta", "Beta-Normal", "Birnbaum-Saunders", "Chisq",
#'     "Exponential", "Exp-Geometric", "Exp-Logarithmic", "Exp-Poisson", "F", "Fisk", "Frechet", "Gamma", "Gaussian", "Gompertz",
#'     "Gumbel", "Kumaraswamy", "Laplace", "Lindley", "Logistic", "Loggamma", "Lognormal", "Lomax", "Makeham", "Maxwell",
#'     "Nakagami", "Paralogistic", "Pareto", "Perks", "Rayleigh", "Rice", "Singh-Maddala", "Skewnormal", "t", "Topp-Leone", "Triangular",
#'     "Uniform", "Weibull".
#'     Please refer to the documentation for each package (either \code{\link[stats]{stats-package}}, \code{\link[VGAM]{VGAM-package}}, or
#'     \code{\link[triangle]{triangle}}) for information on appropriate parameter inputs.
#' @param params a vector of parameters (up to 4) for the desired distribution (keep NULL if \code{fx} supplied instead)
#' @param fx a pdf input as a function of x only, i.e. fx <- function(x) 0.5*(x-1)^2; must return a scalar
#'     (keep NULL if \code{Dist} supplied instead)
#' @param lower the lower support bound for \code{fx}
#' @param upper the upper support bound for \code{fx}
#' @param n the number of random standard normal numbers to use in generating \eqn{y = p(z)} (default = 100)
#' @param seed the seed value for random number generation (default = 1234)
#' @param legend.position the position of the legend
#' @param legend.justification the justification of the legend
#' @param legend.text.size the size of the legend labels
#' @param title.text.size the size of the plot title
#' @param axis.text.size the size of the axes text (tick labels)
#' @param axis.title.size the size of the axes titles
#' @import ggplot2
#' @import stats
#' @import utils
#' @importFrom VGAM dbenini rbenini dbetanorm rbetanorm dbisa rbisa ddagum rdagum dexpgeom rexpgeom dexplog rexplog
#'     dexppois rexppois dfisk rfisk dfrechet rfrechet dgompertz rgompertz dgumbel rgumbel dkumar rkumar dlaplace rlaplace dlind rlind
#'     dlgamma rlgamma dlomax rlomax dmakeham rmakeham dmaxwell rmaxwell dnaka rnaka dparalogistic
#'     rparalogistic dpareto rpareto dperks rperks dgenray rgenray drice rrice dsinmad rsinmad dskewnorm rskewnorm
#'     dtopple rtopple
#' @importFrom triangle dtriangle rtriangle
#' @export
#' @keywords plot, theoretical, pdf, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{calc_theory}},
#'     \code{\link[ggplot2]{ggplot2-package}}, \code{\link[ggplot2]{geom_path}}
#' @return A \code{\link[ggplot2]{ggplot2-package}} object.
#' @references Please see the references for \code{\link[SimMultiCorrData]{plot_cdf}}.
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
                            Dist = c("Benini", "Beta", "Beta-Normal",
                                     "Birnbaum-Saunders", "Chisq", "Dagum",
                                     "Exponential", "Exp-Geometric",
                                     "Exp-Logarithmic", "Exp-Poisson", "F",
                                     "Fisk", "Frechet", "Gamma", "Gaussian",
                                     "Gompertz", "Gumbel", "Kumaraswamy",
                                     "Laplace", "Lindley", "Logistic",
                                     "Loggamma", "Lognormal", "Lomax",
                                     "Makeham", "Maxwell", "Nakagami",
                                     "Paralogistic", "Pareto", "Perks",
                                     "Rayleigh", "Rice", "Singh-Maddala",
                                     "Skewnormal", "t", "Topp-Leone",
                                     "Triangular", "Uniform", "Weibull"),
                            params = NULL, fx = NULL, lower = NULL,
                            upper = NULL, n = 100, seed = 1234,
                            legend.position = c(0.975, 0.9),
                            legend.justification = c(1, 1),
                            legend.text.size = 10, title.text.size = 15,
                            axis.text.size = 10, axis.title.size = 13) {
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
        y[i] <- sigma * (c[1] + c[2] * z[i] + c[3] * z[i]^2 + c[4] * z[i]^3) +
          mu
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
      theme(plot.title = element_text(size = title.text.size, face = "bold",
                                      hjust = 0.5),
            axis.text.x = element_text(size = axis.text.size),
            axis.title.x = element_text(size = axis.title.size),
            axis.text.y = element_text(size = axis.text.size),
            axis.title.y = element_text(size = axis.title.size),
            legend.text = element_text(size = legend.text.size),
            legend.position = legend.position,
            legend.justification = legend.justification) +
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
      D <-
      data.frame(Dist = c("Benini", "Beta", "Beta-Normal", "Birnbaum-Saunders",
                          "Chisq", "Dagum", "Exponential", "Exp-Geometric",
                          "Exp-Logarithmic", "Exp-Poisson", "F", "Fisk",
                          "Frechet", "Gamma", "Gaussian", "Gompertz", "Gumbel",
                          "Kumaraswamy", "Laplace", "Lindley", "Logistic",
                          "Loggamma", "Lognormal", "Lomax", "Makeham",
                          "Maxwell", "Nakagami", "Paralogistic",
                          "Pareto", "Perks", "Rayleigh", "Rice",
                          "Singh-Maddala", "Skewnormal", "t", "Topp-Leone",
                          "Triangular", "Uniform", "Weibull"),
                 pdf = c("dbenini", "dbeta", "dbetanorm", "dbisa", "dchisq",
                         "ddagum", "dexp", "dexpgeom", "dexplog", "dexppois",
                         "df", "dfisk", "dfrechet", "dgamma", "dnorm",
                         "dgompertz", "dgumbel", "dkumar", "dlaplace",
                         "dlind", "dlogis", "dlgamma", "dlnorm",
                         "dlomax", "dmakeham", "dmaxwell", "dnaka",
                         "dparalogistic", "dpareto", "dperks", "dgenray",
                         "drice", "dsinmad", "dskewnorm", "dt", "dtopple",
                         "dtriangle", "dunif", "dweibull"),
                 fx = c("rbenini", "rbeta", "rbetanorm", "rbisa", "rchisq",
                        "rdagum", "rexp", "rexpgeom", "rexplog", "rexppois",
                        "rf", "rfisk", "rfrechet", "rgamma", "rnorm",
                        "rgompertz", "rgumbel", "rkumar", "rlaplace",
                        "rlind", "rlogis", "rlgamma", "rlnorm",
                        "rlomax", "rmakeham", "rmaxwell", "rnaka",
                        "rparalogistic", "rpareto", "rperks", "rgenray",
                        "rrice", "rsinmad", "rskewnorm", "rt", "rtopple",
                        "rtriangle", "runif", "rweibull"),
                 Lower = as.numeric(c(params[1], 0, -Inf, rep(0, 9),
                                      params[1], 0, -Inf, 0, -Inf, 0, -Inf,
                                      0, -Inf, -Inf, rep(0, 6),
                                      params[1], rep(0, 4), -Inf, -Inf, 0,
                                      params[1], params[1], 0)),
                 Upper = as.numeric(c(Inf, 1, rep(Inf, 15), 1, rep(Inf, 17),
                                      1, params[2], params[2], Inf)))
      i <- match(Dist, D$Dist)
      p <- as.character(D$pdf[i])
      for (j in 1:n) {
        if (length(params) == 1) y_fx[j] <- get(p)(x[j], params[1])
        if (length(params) == 2) y_fx[j] <- get(p)(x[j], params[1], params[2])
        if (length(params) == 3) y_fx[j] <- get(p)(x[j], params[1], params[2],
                                                   params[3])
        if (length(params) == 4) y_fx[j] <- get(p)(x[j], params[1], params[2],
                                                   params[3], params[4])
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
      theme(plot.title = element_text(size = title.text.size, face = "bold",
                                      hjust = 0.5),
            axis.text.x = element_text(size = axis.text.size),
            axis.title.x = element_text(size = axis.title.size),
            axis.text.y = element_text(size = axis.text.size),
            axis.title.y = element_text(size = axis.title.size),
            legend.text = element_text(size = legend.text.size),
            legend.position = legend.position,
            legend.justification = legend.justification) +
      scale_colour_manual(name = "", values = c(power_color, target_color),
                          labels = c("Power Method", "Target")) +
      scale_linetype_manual(name ="", values = c(1, target_lty),
                            labels = c("Power Method", "Target"))
    return(plot1)
  }
}
