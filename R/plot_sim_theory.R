#' @title Plot Simulated Data and Target Distribution Data by Name or Function for Continuous or Count Variables
#'
#' @description This plots simulated continuous or count data and overlays data (if \code{overlay} = TRUE) generated from the target
#'     distribution, which is specified by name (plus up to 3 parameters) or pdf function \code{fx} (plus support bounds).
#'     Due to the integration involved in evaluating the cdf using \code{fx}, only continuous \code{fx} may be supplied.  Both are plotted
#'     as histograms.  If a continuous target distribution is specified (\code{cont_var = TRUE}), the simulated data \eqn{y} is
#'     scaled and then transformed (i.e. \eqn{y = sigma * scale(y) + mu}) so that it has the same mean (\eqn{mu}) and variance (\eqn{sigma^2}) as the
#'     target distribution.  If the variable is Negative Binomial, the parameters must be size and success probability (not mu).
#'     It returns a \code{\link[ggplot2]{ggplot2}} object so the user can modify as necessary.
#'     The graph parameters (i.e. \code{title}, \code{power_color}, \code{target_color},
#'     \code{target_lty}) are \code{\link[ggplot2]{ggplot2}} parameters.  It works for valid or invalid power method pdfs.
#' @param sim_y a vector of simulated data
#' @param title the title for the graph (default = "Simulated Data Values")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value)
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value)
#' @param power_color the histogram fill color for the simulated variable (default = "dark blue")
#' @param overlay if TRUE (default), the target distribution is also plotted given either a distribution name (and parameters)
#'     or pdf function fx (with support bounds = lower, upper)
#' @param cont_var TRUE (default) for continuous variables, FALSE for count variables
#' @param target_color the histogram fill color for the target distribution (default = "dark green")
#' @param nbins the number of bins to use when creating the histograms (default = 100)
#' @param Dist name of the distribution. The possible values are: "Beta", "Chisq", "Exponential", "F", "Gamma", "Gaussian",
#'     "Laplace", "Logistic", "Lognormal", "Pareto", "Rayleigh", "t", "Triangular", "Uniform", "Weibull", "Poisson", "Negative_Binomial".
#'     Please refer to the documentation for each package (i.e. \code{\link[stats]{dgamma}})
#'     for information on appropriate parameter inputs.  The pareto (see \code{\link[VGAM]{dpareto}}), generalized
#'     rayleigh (see \code{\link[VGAM]{dgenray}}), and laplace (see \code{\link[VGAM]{dlaplace}}) distributions
#'     come from the \code{\link[VGAM]{VGAM}} package.  The triangular (see \code{\link[triangle]{dtriangle}}) distribution
#'     comes from the \code{\link[triangle]{triangle}} package.
#' @param params a vector of parameters (up to 3) for the desired distribution (keep NULL if \code{fx} supplied instead)
#' @param fx a pdf input as a function of x only, i.e. fx <- function(x) 0.5*(x-1)^2; must return a scalar
#'     (keep NULL if \code{Dist} supplied instead)
#' @param lower the lower support bound for a supplied fx, else keep NULL
#' @param upper the upper support bound for a supplied fx, else keep NULL
#' @param seed the seed value for random number generation (default = 1234)
#' @param sub the number of subdivisions to use in the integration to calculate the cdf from fx; if no result, try increasing
#'     sub (requires longer computation time; default = 1000)
#' @import ggplot2
#' @import stats
#' @import utils
#' @importFrom VGAM dpareto rpareto dgenray rgenray dlaplace rlaplace
#' @importFrom triangle dtriangle rtriangle
#' @export
#' @keywords plot, theoretical, simulated, Fleishman, Headrick
#' @seealso \code{\link[SimMultiCorrData]{calc_theory}},
#'     \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_histogram}}
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
#'                        skews = stcum[3], skurts = stcum[4],
#'                        fifths = stcum[5], sixths = stcum[6],
#'                        n = 10000, seed = seed)
#'
#' # Plot simulated variable (invalid) and data from theoretical distribution
#' plot_sim_theory(sim_y = Logvar1$continuous_variable,
#'                 title = "Invalid Logistic Simulated Data Values",
#'                 overlay = TRUE, Dist = "Logistic", params = c(0, 1),
#'                 seed = seed)
#'
#' # Simulate with the sixth cumulant correction
#' # (valid power method pdf)
#' Logvar2 <- nonnormvar1(method = "Polynomial", means = 0, vars = 1,
#'                        skews = stcum[3], skurts = stcum[4],
#'                        fifths = stcum[5], sixths = stcum[6],
#'                        Six = seq(1.5, 2, 0.05), n = 10000, seed = seed)
#'
#' # Plot simulated variable (valid) and data from theoretical distribution
#' plot_sim_theory(sim_y = Logvar2$continuous_variable,
#'                 title = "Valid Logistic Simulated Data Values",
#'                 overlay = TRUE, Dist = "Logistic", params = c(0, 1),
#'                 seed = seed)
#'
#' # Simulate 2 Negative Binomial distributions and correlation 0.3
#' # using Method 1
#' NBvars <- rcorrvar(k_nb = 2, size = c(10, 15), prob = c(0.4, 0.3),
#'                    rho = matrix(c(1, 0.3, 0.3, 1), 2, 2), seed = seed)
#'
#' # Plot pdfs of 1st simulated variable and theoretical distribution
#' plot_sim_theory(sim_y = NBvars$Neg_Bin_variable[, 1], overlay = TRUE,
#'                 cont_var = FALSE, Dist = "Negative_Binomial",
#'                 params = c(10, 0.4))
#'
#' }
#'
plot_sim_theory <- function(sim_y, title = "Simulated Data Values",
                            ylower = NULL, yupper = NULL,
                            power_color = "dark blue", overlay = TRUE,
                            cont_var = TRUE,
                            target_color = "dark green", nbins = 100,
                            Dist = c("Beta", "Chisq", "Exponential", "F",
                                     "Gamma", "Gaussian", "Laplace",
                                     "Logistic", "Lognormal", "Pareto",
                                     "Rayleigh", "t", "Triangular",
                                     "Uniform", "Weibull", "Poisson",
                                     "Negative_Binomial"),
                            params = NULL, fx = NULL, lower = NULL,
                            upper = NULL, seed = 1234, sub = 1000) {
  if (overlay == FALSE) {
    if (is.null(ylower) & is.null(yupper)) {
      ylower <- min(sim_y)
      yupper <- max(sim_y)
    }
    data <- data.frame(x = 1:length(sim_y), y = sim_y,
                       type = as.factor(rep("sim", length(sim_y))))
    plot1 <- ggplot() + theme_bw() + ggtitle(title) +
      geom_histogram(data = data[data$type == "sim", ],
                     aes_(~y, fill = ~type), bins = nbins) +
      scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 13),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 13),
            legend.text = element_text(size = 10),
            legend.position = c(0.975, 0.9), legend.justification = c(1, 1)) +
      scale_fill_manual(name = "", values = power_color,
                        labels = c("Simulated Variable"))
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
    set.seed(seed)
    n <- length(sim_y)
    if (!is.null(fx)) {
      uni <- runif(n)
      cfx <- function(x, u, FUN = fx) integrate(function(x, FUN = fx) FUN(x),
                                                lower, x, subdivisions = sub,
                                                stop.on.error = FALSE)$value -
        u
      y_fx <- rep(NA, n)
      for (i in 1:n) {
        y_fx[i] <- uniroot(cfx, c(lower, upper), tol = 0.0001, u = uni[i])$root
      }
    }
    if (is.null(fx)) {
      D <- data.frame(Dist = c("Beta", "Chisq", "Exponential", "F", "Gamma",
                               "Gaussian", "Laplace", "Logistic", "Lognormal",
                               "Pareto", "Rayleigh", "t", "Triangular",
                               "Uniform", "Weibull", "Poisson", "Negative_Binomial"),
                      pdf = c("dbeta", "dchisq", "dexp", "df", "dgamma",
                              "dnorm", "dlaplace", "dlogis", "dlnorm",
                              "dpareto", "dgenray", "dt", "dtriangle",
                              "dunif", "dweibull", "dpois", "dnbinom"),
                      fx = c("rbeta", "rchisq", "rexp", "rf", "rgamma",
                             "rnorm", "rlaplace", "rlogis", "rlnorm",
                             "rpareto", "rgenray", "rt", "rtriangle",
                             "runif", "rweibull", "rpois", "rnbinom"),
                      Lower = as.numeric(c(0, 0, 0, 0, 0, -Inf, -Inf, -Inf,
                                           0, params[1], 0, -Inf, params[1],
                                           params[1], 0, 0, 0)),
                      Upper = as.numeric(c(1, Inf, Inf, Inf, Inf, Inf, Inf,
                                           Inf, Inf, Inf, Inf, Inf, params[2],
                                           params[2], Inf, Inf, params[1])))
      i <- match(Dist, D$Dist)
      p <- as.character(D$fx[i])
      x <- 1:n
      y_fx <- numeric(n)
      if (length(params) == 1) y_fx <- get(p)(n, params[1])
      if (length(params) == 2) y_fx <- get(p)(n, params[1], params[2])
      if (length(params) == 3) y_fx <- get(p)(n, params[1], params[2],
                                              params[3])
    }
    data <- data.frame(x = 1:length(sim_y), y = sim_y,
                       type = as.factor(rep("sim", length(sim_y))))
    data2 <- data.frame(x = 1:length(y_fx), y = sort(y_fx),
                        type = as.factor(rep("theory", length(y_fx))))
    data2 <- data.frame(rbind(data, data2))
    plot1 <- ggplot() + theme_bw() + ggtitle(title) +
      geom_histogram(data = data2[data2$type == "sim", ],
                     aes_(~y, fill = ~type), bins = nbins) +
      geom_histogram(data = data2[data2$type == "theory", ],
                     aes_(~y, fill = ~type), bins = nbins) +
      scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10),
            axis.title.x = element_text(size = 13),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 13),
            legend.text = element_text(size = 10),
            legend.position = c(0.975, 0.9), legend.justification = c(1, 1)) +
      scale_fill_manual(name = "", values = c("sim" = power_color,
                                              "theory" = target_color),
                        labels = c("Simulated Variable", "Target Variable"))
    return(plot1)
  }
}
