#' @title Plot Simulated Probability Density Function and Target PDF by Distribution Name or Function for Continuous or Count Variables
#'
#' @description This plots the pdf of simulated continuous or count data and overlays the target pdf (if \code{overlay} = TRUE),
#'     which is specified by distribution name (plus up to 4 parameters) or pdf function \code{fx} (plus support bounds).
#'     If a continuous target distribution is provided (\code{cont_var = TRUE}), the simulated data \eqn{y} is
#'     scaled and then transformed (i.e. \eqn{y = sigma * scale(y) + mu}) so that it has the same mean (\eqn{mu}) and variance (\eqn{sigma^2}) as the
#'     target distribution.  If the variable is Negative Binomial, the parameters must be size and success probability (not mu).
#'     The function returns a \code{\link[ggplot2]{ggplot2-package}} object so the user can modify as necessary.  The graph parameters (i.e. \code{title},
#'     \code{power_color}, \code{target_color}, \code{target_lty}) are \code{\link[ggplot2]{ggplot2-package}} parameters.  It works for valid or invalid power method pdfs.
#' @param sim_y a vector of simulated data
#' @param title the title for the graph (default = "Simulated Probability Density Function")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param power_color the line color for the simulated variable
#' @param overlay if TRUE (default), the target distribution is also plotted given either a distribution name (and parameters)
#'     or pdf function fx (with bounds = ylower, yupper)
#' @param cont_var TRUE (default) for continuous variables, FALSE for count variables
#' @param target_color the line color for the target pdf
#' @param target_lty the line type for the target pdf (default = 2, dashed line)
#' @param Dist name of the distribution. The possible values are: "Benini", "Beta", "Beta-Normal", "Birnbaum-Saunders", "Chisq",
#'     "Exponential", "Exp-Geometric", "Exp-Logarithmic", "Exp-Poisson", "F", "Fisk", "Frechet", "Gamma", "Gaussian", "Gompertz",
#'     "Gumbel", "Kumaraswamy", "Laplace", "Lindley", "Logistic", "Loggamma", "Lognormal", "Lomax", "Makeham", "Maxwell",
#'     "Nakagami", "Paralogistic", "Pareto", "Perks", "Rayleigh", "Rice", "Singh-Maddala", "Skewnormal", "t", "Topp-Leone", "Triangular",
#'     "Uniform", "Weibull", "Poisson", and "Negative_Binomial".
#'     Please refer to the documentation for each package (either \code{\link[stats]{stats-package}}, \code{\link[VGAM]{VGAM-package}}, or
#'     \code{\link[triangle]{triangle}}) for information on appropriate parameter inputs.
#' @param params a vector of parameters (up to 4) for the desired distribution (keep NULL if \code{fx} supplied instead)
#' @param fx a pdf input as a function of x only, i.e. fx <- function(x) 0.5*(x-1)^2; must return a scalar
#'     (keep NULL if \code{Dist} supplied instead)
#' @param lower the lower support bound for \code{fx}
#' @param upper the upper support bound for \code{fx}
#' @param legend.position the position of the legend
#' @param legend.justification the justification of the legend
#' @param legend.text.size the size of the legend labels
#' @param title.text.size the size of the plot title
#' @param axis.text.size the size of the axes text (tick labels)
#' @param axis.title.size the size of the axes titles
#' @import ggplot2
#' @importFrom VGAM dbenini rbenini dbetanorm rbetanorm dbisa rbisa ddagum rdagum dexpgeom rexpgeom dexplog rexplog
#'     dexppois rexppois dfisk rfisk dfrechet rfrechet dgompertz rgompertz dgumbel rgumbel dkumar rkumar dlaplace rlaplace dlind rlind
#'     dlgamma rlgamma dlomax rlomax dmakeham rmakeham dmaxwell rmaxwell dnaka rnaka dparalogistic
#'     rparalogistic dpareto rpareto dperks rperks dgenray rgenray drice rrice dsinmad rsinmad dskewnorm rskewnorm
#'     dtopple rtopple
#' @importFrom triangle dtriangle rtriangle
#' @export
#' @keywords plot, simulated, theoretical, pdf, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{calc_theory}},
#'     \code{\link[ggplot2]{ggplot2-package}}, \code{\link[ggplot2]{geom_path}}, \code{\link[ggplot2]{geom_density}}
#' @return A \code{\link[ggplot2]{ggplot2-package}} object.
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
#'                        skews = stcum[3], skurts = stcum[4],
#'                        fifths = stcum[5], sixths = stcum[6],
#'                        n = 10000, seed = seed)
#'
#' # Plot pdfs of simulated variable (invalid) and theoretical distribution
#' plot_sim_pdf_theory(sim_y = Logvar1$continuous_variable,
#'                     title = "Invalid Logistic Simulated PDF",
#'                     overlay = TRUE, Dist = "Logistic", params = c(0, 1))
#'
#' # Simulate with the sixth cumulant correction
#' # (valid power method pdf)
#' Logvar2 <- nonnormvar1(method = "Polynomial", means = 0, vars = 1,
#'                        skews = stcum[3], skurts = stcum[4],
#'                        fifths = stcum[5], sixths = stcum[6],
#'                        Six = seq(1.5, 2, 0.05), n = 10000, seed = seed)
#'
#' # Plot pdfs of simulated variable (invalid) and theoretical distribution
#' plot_sim_pdf_theory(sim_y = Logvar2$continuous_variable,
#'                     title = "Valid Logistic Simulated PDF",
#'                     overlay = TRUE, Dist = "Logistic", params = c(0, 1))
#'
#' # Simulate 2 Negative Binomial distributions and correlation 0.3
#' # using Method 1
#' NBvars <- rcorrvar(k_nb = 2, size = c(10, 15), prob = c(0.4, 0.3),
#'                   rho = matrix(c(1, 0.3, 0.3, 1), 2, 2), seed = seed)
#'
#' # Plot pdfs of 1st simulated variable and theoretical distribution
#' plot_sim_pdf_theory(sim_y = NBvars$Neg_Bin_variable[, 1], overlay = TRUE,
#'                     cont_var = FALSE, Dist = "Negative_Binomial",
#'                     params = c(10, 0.4))
#'
#' }
#'
plot_sim_pdf_theory <- function(sim_y,
                       title = "Simulated Probability Density Function",
                       ylower = NULL, yupper = NULL, power_color = "dark blue",
                       overlay = TRUE, cont_var = TRUE,
                       target_color = "dark green", target_lty = 2,
                       Dist = c("Benini", "Beta", "Beta-Normal",
                                "Birnbaum-Saunders", "Chisq", "Dagum",
                                "Exponential", "Exp-Geometric",
                                "Exp-Logarithmic", "Exp-Poisson", "F",
                                "Fisk", "Frechet", "Gamma", "Gaussian",
                                "Gompertz", "Gumbel", "Kumaraswamy",
                                "Laplace", "Lindley", "Logistic",
                                "Loggamma", "Lognormal", "Lomax", "Makeham",
                                "Maxwell", "Nakagami", "Paralogistic",
                                "Pareto", "Perks", "Rayleigh", "Rice",
                                "Singh-Maddala", "Skewnormal", "t",
                                "Topp-Leone", "Triangular", "Uniform",
                                "Weibull", "Poisson", "Negative_Binomial"),
                       params = NULL, fx = NULL, lower = NULL, upper = NULL,
                       legend.position = c(0.975, 0.9),
                       legend.justification = c(1, 1), legend.text.size = 10,
                       title.text.size = 15, axis.text.size = 10,
                       axis.title.size = 13) {
  if (overlay == FALSE) {
    if (is.null(ylower) & is.null(yupper)) {
      ylower <- min(sim_y)
      yupper <- max(sim_y)
    }
    data <- data.frame(x = 1:length(sim_y), y = sim_y,
                       type = as.factor(rep("sim", length(sim_y))))
    plot1 <- ggplot() + theme_bw() + ggtitle(title) +
      geom_density(data = data, aes_(x = ~y, colour = "Density",
                                     lty = ~type)) +
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
      scale_linetype_manual(name = "", values = 1,
                            labels = "Simulated Variable") +
      scale_colour_manual(name = "", values = c(power_color),
                          labels = "Simulated Variable")
    return(plot1)
  }
  if (overlay == TRUE) {
    if (cont_var == TRUE) {
      if (!is.null(fx)) {
        theory <- calc_theory(fx = fx, lower = lower, upper = upper)
      }
      if (is.null(fx)) {
        theory <- calc_theory(Dist = Dist, params = params)
      }
      mu <- theory[1]
      sigma <- theory[2]
      sim_y <- sigma * scale(sim_y) + mu
    }
    if (is.null(ylower) & is.null(yupper)) {
      ylower <- min(sim_y)
      yupper <- max(sim_y)
    }
    x <- sim_y
    y_fx <- numeric(length(x))
    if (!is.null(fx)) {
      for (j in 1:length(x)) {
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
                          "Loggamma", "Lognormal", "Lomax",
                          "Makeham", "Maxwell", "Nakagami", "Paralogistic",
                          "Pareto", "Perks", "Rayleigh", "Rice",
                          "Singh-Maddala", "Skewnormal", "t", "Topp-Leone",
                          "Triangular", "Uniform", "Weibull", "Poisson",
                          "Negative_Binomial"),
                 pdf = c("dbenini", "dbeta", "dbetanorm", "dbisa", "dchisq",
                         "ddagum", "dexp", "dexpgeom", "dexplog", "dexppois",
                         "df", "dfisk", "dfrechet", "dgamma", "dnorm",
                         "dgompertz", "dgumbel", "dkumar", "dlaplace",
                         "dlind", "dlogis", "dlgamma", "dlnorm",
                         "dlomax", "dmakeham", "dmaxwell", "dnaka",
                         "dparalogistic", "dpareto", "dperks", "dgenray",
                         "drice", "dsinmad", "dskewnorm", "dt", "dtopple",
                         "dtriangle", "dunif", "dweibull", "dpois", "dnbinom"),
                 fx = c("rbenini", "rbeta", "rbetanorm", "rbisa", "rchisq",
                        "rdagum", "rexp", "rexpgeom", "rexplog", "rexppois",
                        "rf", "rfisk", "rfrechet", "rgamma", "rnorm",
                        "rgompertz", "rgumbel", "rkumar", "rlaplace",
                        "rlind", "rlogis", "rlgamma", "rlnorm",
                        "rlomax", "rmakeham", "rmaxwell", "rnaka",
                        "rparalogistic", "rpareto", "rperks", "rgenray",
                        "rrice", "rsinmad", "rskewnorm", "rt", "rtopple",
                        "rtriangle", "runif", "rweibull", "rpois", "rnbinom"),
                 Lower = as.numeric(c(params[1], 0, -Inf, rep(0, 9),
                                      params[1], 0, -Inf, 0, -Inf, 0, -Inf,
                                      0, -Inf, -Inf, rep(0, 6),
                                      params[1], rep(0, 4), -Inf, -Inf, 0,
                                      params[1], params[1], 0, 0, 0)),
                 Upper = as.numeric(c(Inf, 1, rep(Inf, 15), 1, rep(Inf, 17),
                                      1, params[2], params[2], Inf, Inf, Inf)))
      i <- match(Dist, D$Dist)
      p <- as.character(D$pdf[i])
      for (j in 1:length(x)) {
        if (length(params) == 1) y_fx[j] <- get(p)(x[j], params[1])
        if (length(params) == 2) y_fx[j] <- get(p)(x[j], params[1], params[2])
        if (length(params) == 3) y_fx[j] <- get(p)(x[j], params[1], params[2],
                                                   params[3])
        if (length(params) == 4) y_fx[j] <- get(p)(x[j], params[1], params[2],
                                                   params[3], params[4])
      }
    }
    data <- data.frame(x = 1:length(sim_y), y = sim_y,
                       type = as.factor(rep("sim", length(sim_y))))
    data2 <- data.frame(x = x, y = y_fx,
                        type = as.factor(rep("theory", length(y_fx))))
    data2 <- data.frame(rbind(data, data2))
    plot1 <- ggplot() + theme_bw() + ggtitle(title) +
      geom_density(data = data2[data2$type == "sim", ],
                   aes_(x = ~y, colour = ~type, lty = ~type)) +
      geom_line(data = data2[data2$type == "theory", ],
                aes_(x = ~x, y = ~y, colour = ~type, lty = ~type)) +
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
      scale_linetype_manual(name = "", values = c(1, target_lty),
                            labels = c("Simulated Variable", "Target")) +
      scale_colour_manual(name = "", values = c(power_color, target_color),
                          labels = c("Simulated Variable", "Target"))
    return(plot1)
  }
}
