#' @title Generation of Correlated Ordinal, Continuous, Poisson, and/or Negative Binomial Variables: Method 1
#'
#' @description This function simulates \code{k_cat} ordinal, \code{k_cont} continuous, \code{k_pois} Poisson, and/or \code{k_nb}
#'     Negative Binomial variables with
#'     a specified correlation matrix \code{rho}.  The variables are generated from multivariate normal variables with intermediate correlation
#'     matrix \code{Sigma}, calculated by \code{\link[SimMultiCorrData]{findintercorr}}, and then transformed.  The \emph{ordering} of the
#'     variables in \code{rho} must be \emph{ordinal} (r >= 2 categories), \emph{continuous}, \emph{Poisson}, and \emph{Negative Binomial}
#'     (note that it is possible for \code{k_cat}, \code{k_cont}, \code{k_pois}, and/or \code{k_nb} to be 0).  The vignette
#'     \bold{Overall Workflow for Data Simulation} provides a detailed example discussing the step-by-step simulation process and comparing
#'     methods 1 and 2.
#'
#' @section Overview of Method 1:
#'     The intermediate correlations used in method 1 are more simulation based than those in method 2, which means that accuracy
#'     increases with sample size and the number of repetitions.  In addition, specifying the seed allows for reproducibility.  In
#'     addition, method 1 differs from method 2 in the following ways:
#'
#'     1) The intermediate correlation for \bold{count variables} is based on the method of Yahav & Shmueli (2012), which uses a
#'     simulation based, logarithmic transformation of the target correlation.  This method becomes less accurate as the variable mean
#'     gets closer to zero.
#'
#'     2) The \bold{ordinal - count variable} correlations are based on an extension of the method of Amatya & Demirtas (2015), in which
#'     the correlation correction factor is the product of the upper Frechet-Hoeffding bound on the correlation between the count
#'     variable and the normal variable used to generate it and a simulated upper bound on the correlation between an ordinal variable
#'     and the normal variable used to generate it (see Demirtas & Hedeker, 2011).
#'
#'     3) The \bold{continuous - count variable} correlations are based on an extension of the methods of Amatya & Demirtas (2015) and
#'     Demirtas et al. (2012), in which the correlation correction factor is the product of the upper Frechet-Hoeffding bound
#'     on the correlation between the count variable and the normal variable used to generate it and the power method correlation
#'     between the continuous variable and the normal variable used to generate it (see Headrick & Kowalchuk, 2007).  The
#'     intermediate correlations are the ratio of the target correlations to the correction factor.
#'
#'     Please see the \bold{Comparison of Method 1 and Method 2} vignette for more information and an step-by-step overview of the
#'     simulation process.
#'
#' @section Reasons for Function Errors:
#'     The most likely cause for function errors is that no solutions to \code{\link[SimMultiCorrData]{fleish}} or
#'     \code{\link[SimMultiCorrData]{poly}} converged when using \code{\link[SimMultiCorrData]{find_constants}}.  If this happens,
#'     the simulation will stop.  It may help to first use \code{\link[SimMultiCorrData]{find_constants}} for each continuous variable to
#'     determine if a vector of sixth cumulant correction values is needed.  The solutions can be used as starting values (see \code{cstart} below).
#'     In addition, the kurtosis may be outside the region of possible values.  There is an associated lower boundary for kurtosis associated
#'     with a given skew (for Fleishman's method) or skew and fifth and sixth cumulants (for Headrick's method).  Use
#'     \code{\link[SimMultiCorrData]{calc_lower_skurt}} to determine the boundary for a given set of cumulants.
#'
#'     In addition, as mentioned above, the feasibility of the final correlation matrix rho, given the
#'     distribution parameters, should be checked first using \code{\link[SimMultiCorrData]{valid_corr}}.  This function either checks
#'     if a given \code{rho} is plausible or returns the lower and upper final correlation limits.  It should be noted that even if a target
#'     correlation matrix is within the "plausible range," it still may not be possible to achieve the desired matrix.  This happens most
#'     frequently when generating ordinal variables (r >= 2 categories).  The error loop frequently fixes these problems.
#'
#' @param n the sample size (i.e. the length of each simulated variable; default = 10000)
#' @param k_cont the number of continuous variables (default = 0)
#' @param k_cat the number of ordinal (r >= 2 categories) variables (default = 0)
#' @param k_pois the number of Poisson variables (default = 0)
#' @param k_nb the number of Negative Binomial variables (default = 0)
#' @param method the method used to generate the k_cont continuous variables.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param means a vector of means for the k_cont continuous variables (i.e. = rep(0, k_cont))
#' @param vars a vector of variances (i.e. = rep(1, k_cont))
#' @param skews a vector of skewness values (i.e. = rep(0, k_cont))
#' @param skurts a vector of standardized kurtoses (kurtosis - 3, so that normal variables have a value of 0; i.e. = rep(0, k_cont))
#' @param fifths a vector of standardized fifth cumulants (not necessary for \code{method} = "Fleishman"; i.e. = rep(0, k_cont))
#' @param sixths a vector of standardized sixth cumulants (not necessary for \code{method} = "Fleishman"; i.e. = rep(0, k_cont))
#' @param Six a list of vectors of correction values to add to the sixth cumulants if no valid pdf constants are found,
#'     ex: \code{Six = list(seq(0.01, 2,by = 0.01), seq(1, 10,by = 0.5))}; if no correction is desired for variable Y_i, set set the i-th list
#'     component equal to NULL
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1; default = list());
#'     for binary variables, these should be input the same as for ordinal variables with more than 2 categories (i.e. the user-specified
#'     probability is the probability of the 1st category, which has the smaller support value)
#' @param support a list of length equal to \code{k_cat}; the i-th element is a vector containing the r ordered support values;
#'     if not provided (i.e. \code{support = list()}), the default is for the i-th element to be the
#'     vector 1, ..., r
#' @param nrand the number of random numbers to generate in calculating intermediate correlations (default = 10000)
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{dpois}})
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param Sigma an intermediate correlation matrix to use if the user wants to provide one (default = NULL)
#' @param rho the target correlation matrix (\emph{must be ordered ordinal, continuous, Poisson, Negative Binomial}; default = NULL)
#' @param cstart a list containing initial values for root-solving algorithm used in \code{\link[SimMultiCorrData]{find_constants}}
#'     (see \code{\link[BB]{multiStart}} for \code{method} = "Fleishman" or \code{\link[nleqslv]{nleqslv}} for \code{method} = "Polynomial").
#'     If user specified, each list element must be input as a matrix. If no starting values are specified for a given continuous
#'     variable, that list element should be NULL.  If NULL and all 4 standardized cumulants (rounded to 3 digits) are within
#'     0.01 of those in Headrick's common distribution table (see \code{\link[SimMultiCorrData]{Headrick.dist}}
#'     data), uses his constants as starting values; else, generates n sets of random starting values from
#'     uniform distributions.
#' @param seed the seed value for random number generation (default = 1234)
#' @param errorloop if TRUE, uses \code{\link[SimMultiCorrData]{error_loop}} to attempt to correct the final correlation
#'     (default = FALSE)
#' @param epsilon the maximum acceptable error between the final and target correlation matrices (default = 0.001)
#'     in the calculation of ordinal intermediate correlations with \code{\link[SimMultiCorrData]{ordnorm}} or in the error loop
#' @param maxit the maximum number of iterations to use (default = 1000) in the calculation of ordinal
#'     intermediate correlations with \code{\link[SimMultiCorrData]{ordnorm}} or in the error loop
#' @param extra_correct if TRUE, within each variable pair, if the maximum correlation error is still greater than 0.1, the intermediate
#'     correlation is set equal to the target correlation (with the assumption that the calculated final correlation will be
#'     less than 0.1 away from the target)
#' @importFrom psych describe
#' @import stats
#' @import utils
#' @importFrom Matrix nearPD
#' @import GenOrd
#' @import BB
#' @import nleqslv
#' @export
#' @keywords simulation, continuous, ordinal, Poisson, Negative Binomial, Fleishman, Headrick, method1
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{findintercorr}},
#'     \code{\link[BB]{multiStart}}, \code{\link[nleqslv]{nleqslv}}
#' @return A list whose components vary based on the type of simulated variables.  Simulated variables are returned as data.frames:
#' @return If \bold{ordinal variables} are produced:
#'
#'     \code{ordinal_variables} the generated ordinal variables,
#'
#'     \code{summary_ordinal} a list, where the i-th element contains a data.frame with column 1 = target cumulative probabilities and
#'     column 2 = simulated cumulative probabilities for ordinal variable Y_i
#' @return If \bold{continuous variables} are produced:
#'
#'     \code{constants} a data.frame of the constants,
#'
#'     \code{continuous_variables} the generated continuous variables,
#'
#'     \code{summary_continuous} a data.frame containing a summary of each variable,
#'
#'     \code{summary_targetcont} a data.frame containing a summary of the target variables,
#'
#'     \code{sixth_correction} a vector of sixth cumulant correction values,
#'
#'     \code{valid.pdf} a vector where the i-th element is "TRUE" if the constants for the i-th continuous variable generate a valid pdf,
#'                 else "FALSE"
#' @return If \bold{Poisson variables} are produced:
#'
#'     \code{Poisson_variables} the generated Poisson variables,
#'
#'     \code{summary_Poisson} a data.frame containing a summary of each variable
#' @return If \bold{Negative Binomial variables} are produced:
#'
#'     \code{Neg_Bin_variables} the generated Negative Binomial variables,
#'
#'     \code{summary_Neg_Bin} a data.frame containing a summary of each variable
#' @return Additionally, the following elements:
#'
#'     \code{correlations} the final correlation matrix,
#'
#'     \code{Sigma1} the intermediate correlation before the error loop,
#'
#'     \code{Sigma2} the intermediate correlation matrix after the error loop,
#'
#'     \code{Constants_Time} the time in minutes required to calculate the constants,
#'
#'     \code{Intercorrelation_Time} the time in minutes required to calculate the intermediate correlation matrix,
#'
#'     \code{Error_Loop_Time} the time in minutes required to use the error loop,
#'
#'     \code{Simulation_Time} the total simulation time in minutes,
#'
#'     \code{niter} a matrix of the number of iterations used for each variable in the error loop,
#'
#'     \code{maxerr} the maximum final correlation error (from the target rho).
#'
#'     If a particular element is not required, the result is NULL for that element.
#' @references Ferrari PA, Barbiero A (2012). Simulating ordinal data, Multivariate Behavioral Research, 47(4): 566-589.
#'
#' Barbiero A, Ferrari PA (2015). GenOrd: Simulation of Discrete Random Variables with Given
#' Correlation Matrix and Marginal Distributions. R package version 1.4.0. \url{https://CRAN.R-project.org/package=GenOrd}
#'
#' Higham N (2002). Computing the nearest correlation matrix - a problem from finance; IMA Journal of Numerical Analysis 22: 329-343.
#'
#' Olsson U, Drasgow F, & Dorans NJ (1982). The Polyserial Correlation Coefficient. Psychometrika, 47(3): 337-47.
#'     \doi{10.1007/BF02294164}.
#'
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis 40(4):685-711
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43, 521-532.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64, 25-35.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77, 229-249.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3, 65-71.
#'
#' Demirtas H, Hedeker D, & Mermelstein RJ (2012). Simulation of massive public health data by power polynomials.
#'     Statistics in Medicine 31:27, 3337-3346.
#'
#' Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15): 3129-39.
#'
#' Amatya A & Demirtas H (2016). PoisNor: Simultaneous Generation of Multivariate Data with Poisson and Normal Marginals.
#'     R package version 1.1. \url{https://CRAN.R-project.org/package=PoisNor}
#'
#' Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1): 91-102. \doi{10.1002/asmb.901}.
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2): 104-109.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Frechet M.  Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA.  1951;14:53-77.
#'
#' Varadhan R, Gilbert PD (2009). BB: An R Package for Solving a Large System of Nonlinear Equations and for
#'     Optimizing a High-Dimensional Nonlinear Objective Function, J. Statistical Software, 32:4,
#'     \url{http://www.jstatsoft.org/v32/i04/}
#'
#' Berend Hasselman (2017). nleqslv: Solve Systems of Nonlinear Equations. R package version 3.2.
#'     \url{https://CRAN.R-project.org/package=nleqslv}
#' @examples \dontrun{
#'
#' # Binary, Ordinal, Continuous, Poisson, and Negative Binomial Variables
#'
#' options(scipen = 999)
#' seed <- 1234
#' n <- 10000
#'
#' # Continuous Distributions: Normal, t (df = 10), Chisq (df = 4),
#' #                           Beta (a = 4, b = 2), Gamma (a = 4, b = 4)
#' Dist <- c("Gaussian", "t", "Chisq", "Beta", "Gamma")
#'
#' # calculate standardized cumulants
#'
#' M1 <- calc_theory(Dist = "Gaussian", params = c(0, 1))
#' M2 <- calc_theory(Dist = "t", params = 10)
#' M3 <- calc_theory(Dist = "Chisq", params = 4)
#' M4 <- calc_theory(Dist = "Beta", params = c(4, 2))
#' M5 <- calc_theory(Dist = "Gamma", params = c(4, 4))
#' M <- cbind(M1, M2, M3, M4, M5)
#' M <- round(M[-c(1:2),], digits = 6)
#' colnames(M) <- Dist
#' rownames(M) <- c("skew", "skurtosis", "fifth", "sixth")
#' means <- rep(0, length(Dist))
#' vars <- rep(1, length(Dist))
#'
#' # Binary and Ordinal Distributions
#' marginal <- list(0.3, 0.4, c(0.1, 0.5), c(0.3, 0.6, 0.9),
#'                  c(0.2, 0.4, 0.7, 0.8))
#' support <- list()
#'
#' # Poisson Distributions
#' lam <- c(1, 5, 10)
#'
#' # Negative Binomial Distributions
#' size <- c(3, 6)
#' prob <- c(0.2, 0.8)
#'
#' ncat <- length(marginal)
#' ncont <- ncol(M)
#' npois <- length(lam)
#' nnb <- length(size)
#'
#' # Create correlation matrix from a uniform distribution (-0.8, 0.8)
#' set.seed(seed)
#' Rey <- diag(1, nrow = (ncat + ncont + npois + nnb))
#' for (i in 1:nrow(Rey)) {
#'   for (j in 1:ncol(Rey)) {
#'     if (i > j) Rey[i, j] <- runif(1, -0.8, 0.8)
#'     Rey[j, i] <- Rey[i, j]
#'   }
#' }
#'
#' # Test for positive-definiteness
#' library(Matrix)
#' if(min(eigen(Rey, symmetric = TRUE)$values) < 0) {
#'   Rey <- as.matrix(nearPD(Rey, corr = T, keepDiag = T)$mat)
#' }
#'
#' # Make sure Rey is within upper and lower correlation limits
#' valid <- valid_corr(k_cat = ncat, k_cont = ncont, k_pois = npois,
#'                     k_nb = nnb, method = "Polynomial", means = means,
#'                     vars = vars, skews = M[1, ], skurts = M[2, ],
#'                     fifths = M[3, ], sixths = M[4, ], marginal = marginal,
#'                     lam = lam, size = size, prob = prob, rho = Rey,
#'                     seed = seed)
#'
#' # Simulate variables without error loop
#' A <- rcorrvar(n = 10000, k_cont = ncont, k_cat = ncat, k_pois = npois,
#'               k_nb = nnb, method = "Polynomial", means = means, vars = vars,
#'               skews = M[1, ], skurts = M[2, ], fifths = M[3, ],
#'               sixths = M[4, ], marginal = marginal, lam = lam, size = size,
#'               prob = prob, rho = Rey, seed = seed)
#'
#' # Look at the maximum correlation error
#' A$maxerr
#' Acorr_error = round(A$correlations - Rey, 6)
#'
#' # interquartile-range of correlation errors
#' quantile(as.numeric(Acorr_error), 0.25)
#' quantile(as.numeric(Acorr_error), 0.75)
#'
#' # Simulate variables with error loop (using default settings of
#' # epsilon = 0.001 and maxit = 1000)
#' B <- rcorrvar(n = 10000, k_cont = ncont, k_cat = ncat, k_pois = npois,
#'               k_nb = nnb, method = "Polynomial", means = means, vars = vars,
#'               skews = M[1, ], skurts = M[2, ], fifths = M[3, ],
#'               sixths = M[4, ], marginal = marginal, lam = lam, size = size,
#'               prob = prob, rho = Rey, seed = seed, errorloop = TRUE)
#'
#' # Look at the maximum correlation error
#' B$maxerr
#' Bcorr_error = round(B$correlations - Rey, 6)
#'
#' # interquartile-range of correlation errors
#' quantile(as.numeric(Bcorr_error), 0.25)
#' quantile(as.numeric(Bcorr_error), 0.75)
#'
#' # Look at results
#' # Ordinal variables
#' B$summary_ordinal
#'
#' # Continuous variables
#' round(B$constants, 6)
#' round(B$summary_continuous, 6)
#' round(B$summary_targetcont, 6)
#' B$valid.pdf
#'
#' # Count variables
#' B$summary_Poisson
#' B$summary_Neg_Bin
#'
#' # Generate Plots
#'
#' # t (df = 10) (2nd continuous variable)
#' # 1) Simulated Data CDF (find cumulative probability up to y = 0.5)
#' plot_sim_cdf(B$continuous_variables[, 2], calc_cprob = TRUE, delta = 0.5)
#'
#' # 2) Simulated Data and Target Distribution PDFs
#' plot_sim_pdf_theory(B$continuous_variables[, 2], Dist = "t", params = 10)
#'
#' # 3) Simulated Data and Target Distribution
#' plot_sim_theory(B$continuous_variables[, 2], Dist = "t", params = 10)
#'
#' # Chisq (df = 4) (3rd continuous variable)
#' # 1) Simulated Data CDF (find cumulative probability up to y = 0.5)
#' plot_sim_cdf(B$continuous_variables[, 3], calc_cprob = TRUE, delta = 0.5)
#'
#' # 2) Simulated Data and Target Distribution PDFs
#' plot_sim_pdf_theory(B$continuous_variables[, 3], Dist = "Chisq", params = 4)
#'
#' # 3) Simulated Data and Target Distribution
#' plot_sim_theory(B$continuous_variables[, 3], Dist = "Chisq", params = 4)
#'
#' }
rcorrvar <- function(n = 10000, k_cont = 0, k_cat = 0, k_pois = 0, k_nb = 0,
                     method = c("Fleishman", "Polynomial"),
                     means =  NULL, vars =  NULL, skews =  NULL,
                     skurts =  NULL, fifths =  NULL, sixths =  NULL,
                     Six = list(), marginal = list(), support = list(),
                     nrand = 100000, lam  =  NULL,
                     size = NULL, prob = NULL, mu = NULL, Sigma = NULL,
                     rho = NULL, cstart = NULL, seed = 1234,
                     errorloop = FALSE, epsilon = 0.001,  maxit = 1000,
                     extra_correct = TRUE) {
  start.time <- Sys.time()
  k <- k_cat + k_cont + k_pois + k_nb
  if (k_cat > 0 & k_cat != length(marginal))
    stop("Length of marginal does not match the number of Ordinal variables!")
  if (k_cont > 0 & k_cont != length(means))
    stop("Length of means does not match the number of Continuous variables!")
  if (k_cont > 0 & k_cont != length(vars))
    stop("Length of variances does not match the number of Continuous
         variables!")
  if (k_pois > 0 & k_pois != length(lam))
    stop("Length of lam does not match the number of Poisson variables!")
  if (k_pois > 0 & sum(lam < 0) > 0)
    stop("Lambda values cannnot be negative!")
  if (k_nb > 0 & k_nb != length(size))
    stop("Length of size does not match the number of Negative Binomial
         variables!")
  if (k_nb > 0 & length(prob) > 1 & length(mu) > 1)
    stop("Either give success probabilities or means for Negative Binomial
         variables!")
  SixCorr <- numeric(k_cont)
  Valid.PDF <- numeric(k_cont)
  if (ncol(rho) != k)
    stop("Correlation matrix dimension does not match the number of
         variables!")
  if (!isSymmetric(rho) | min(eigen(rho, symmetric = TRUE)$values) < 0 |
      !all(diag(rho) == 1)) stop("Correlation matrix not valid!")
  start.time.constants <- Sys.time()
  if (k_cont >= 1) {
    if (method == "Fleishman") {
      constants <- matrix(1, nrow = k_cont, ncol = 4)
      colnames(constants) <- c("c0", "c1", "c2", "c3")
    }
    if (method == "Polynomial") {
      constants <- matrix(1, nrow = k_cont, ncol = 6)
      colnames(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5")
    }
    for (i in 1:k_cont) {
      if (length(Six) == 0) {
        if (length(cstart) == 0) {
          cons <-
            suppressWarnings(find_constants(method = method, skews = skews[i],
                                            skurts = skurts[i],
                                            fifths = fifths[i],
                                            sixths = sixths[i], Six = Six,
                                            cstart = cstart, n = 25,
                                            seed = seed))
        } else {
          cons <-
            suppressWarnings(find_constants(method = method, skews = skews[i],
                                            skurts = skurts[i],
                                            fifths = fifths[i],
                                            sixths = sixths[i], Six = Six,
                                            cstart = cstart[[i]], n = 25,
                                            seed = seed))
        }
      } else {
        if (length(cstart) == 0) {
          cons <-
            suppressWarnings(find_constants(method = method, skews = skews[i],
                                            skurts = skurts[i],
                                            fifths = fifths[i],
                                            sixths = sixths[i], Six = Six[[i]],
                                            cstart = cstart, n = 25,
                                            seed = seed))
        } else {
          cons <-
            suppressWarnings(find_constants(method = method, skews = skews[i],
                                            skurts = skurts[i],
                                            fifths = fifths[i],
                                            sixths = sixths[i], Six = Six[[i]],
                                            cstart = cstart[[i]], n = 25,
                                            seed = seed))
        }
      }
      if (length(cons) == 1) {
        stop(paste("Constants can not be found for continuous variable ", i,
                   ".", sep = ""))
      }
      con_solution <- cons$constants
      SixCorr[i] <- ifelse(is.null(cons$SixCorr1), NA, cons$SixCorr1)
      Valid.PDF[i] <- cons$valid
      constants[i, ] <- con_solution
      cat("\n", "Constants: Distribution ", i, " \n")
    }
  }
  stop.time.constants <- Sys.time()
  if (k_cat > 0) {
    if (!all(unlist(lapply(marginal, function(x) (sort(x) == x & min(x) > 0 &
                                                  max(x) < 1)))))
      stop("Error in given marginal distributions!")
    len <- length(support)
    kj <- numeric(k_cat)
    for (i in 1:k_cat) {
      kj[i] <- length(marginal[[i]]) + 1
      if (len == 0) {
        support[[i]] <- 1:kj[i]
      }
    }
  }
  start.time.intercorr <- Sys.time()
  if (is.null(Sigma)) {
    Sigma <- findintercorr(n = n, k_cont = k_cont, k_cat = k_cat,
                           k_pois = k_pois, k_nb = k_nb, method = method,
                           constants = constants, marginal = marginal,
                           support = support, nrand = nrand, lam = lam,
                           size = size, prob = prob, mu = mu, rho = rho,
                           seed = seed, epsilon = epsilon, maxit = maxit)
  }
  if (min(eigen(Sigma, symmetric = TRUE)$values) < 0) {
    warning("Intermediate correlation matrix is not positive definite.
            Nearest positive definite matrix is used!")
    Sigma <- as.matrix(nearPD(Sigma, corr = T, keepDiag = T)$mat)
  }
  if (!isSymmetric(Sigma) | min(eigen(Sigma, symmetric = TRUE)$values) < 0 |
      !all(diag(Sigma) == 1))
    stop("Calculated intermediate correlation matrix not valid!")
  stop.time.intercorr <- Sys.time()
  eig <- eigen(Sigma, symmetric = TRUE)
  sqrteigval <- diag(sqrt(eig$values), nrow = nrow(Sigma), ncol = ncol(Sigma))
  eigvec <- eig$vectors
  fry <- eigvec %*% sqrteigval
  set.seed(seed)
  X <- matrix(rnorm(k * n), n)
  X <- scale(X, TRUE, FALSE)
  X <- X %*% svd(X, nu = 0)$v
  X <- scale(X, FALSE, TRUE)
  X <- fry %*% t(X)
  X <- t(X)
  if (k_cat > 0) {
    X_cat <- matrix(X[, 1:k_cat], nrow = n, ncol = k_cat, byrow = F)
  }
  if (k_cont > 0) {
    X_cont <- matrix(X[, (k_cat + 1):(k_cat + k_cont)], nrow = n,
                     ncol = k_cont, byrow = F)
  }
  if (k_pois > 0) {
    X_pois <- matrix(X[, (k_cat + k_cont + 1):(k_cat + k_cont + k_pois)],
                     nrow = n, ncol = k_pois, byrow = F)
  }
  if (k_nb > 0) {
    X_nb <- matrix(X[, (k_cat + k_cont + k_pois + 1):(k_cat + k_cont +
                                                        k_pois + k_nb)],
                   nrow = n, ncol = k_nb, byrow = F)
  }
  Y_cat <- NA
  Y <- NA
  Yb <- NA
  Y_pois <- NA
  Y_nb <- NA
  if (k_cat > 0) {
    Y_cat <- matrix(1, nrow = n, ncol = k_cat)
    for (i in 1:length(marginal)) {
      Y_cat[, i] <- as.integer(cut(X_cat[, i], breaks = c(min(X_cat[, i]) - 1,
                               qnorm(marginal[[i]]), max(X_cat[, i])  +  1)))
      Y_cat[, i] <- support[[i]][Y_cat[, i]]
    }
  }
  if (k_cont > 0) {
    Y <- matrix(1, nrow = n, ncol = k_cont)
    Yb <- matrix(1, nrow = n, ncol = k_cont)
    for (i in 1:k_cont) {
      if (method == "Fleishman") {
        Y[, i] <- constants[i, 1] + constants[i, 2] * X_cont[, i] +
          constants[i, 3] * X_cont[, i]^2 + constants[i, 4] * X_cont[, i]^3
      }
      if (method == "Polynomial") {
        Y[, i] <- constants[i, 1] + constants[i, 2] * X_cont[, i] +
          constants[i, 3] * X_cont[, i]^2 + constants[i, 4] * X_cont[, i]^3 +
          constants[i, 5] * X_cont[, i]^4 + constants[i, 6] * X_cont[, i]^5
      }
      Yb[, i] <- means[i] + sqrt(vars[i]) * Y[, i]
    }
  }
  if (k_pois > 0) {
    Y_pois <- matrix(1, nrow = n, ncol = k_pois)
    for (i in 1:k_pois) {
      Y_pois[, i] <- qpois(pnorm(X_pois[, i]), lam[i])
    }
  }
  if (k_nb > 0) {
    Y_nb <- matrix(1, nrow = n, ncol = k_nb)
    if (length(prob) > 0) {
      for (i in 1:k_nb) {
        Y_nb[, i] <- qnbinom(pnorm(X_nb[, i]), size[i], prob[i])
      }
    }
    if (length(mu) > 0) {
      for (i in 1:k_nb) {
        Y_nb[, i] <- qnbinom(pnorm(X_nb[, i]), size[i], mu = mu[i])
      }
    }
  }
  rho_calc <- calc_final_corr(k_cat = k_cat, k_cont = k_cont, k_pois = k_pois,
                              k_nb = k_nb, Y_cat = Y_cat, Yb = Yb,
                              Y_pois = Y_pois, Y_nb = Y_nb)
  Sigma1 <- Sigma
  start.time.error <- Sys.time()
  rho0 <- rho
  niter <- diag(0, k, k)
  emax <- max(abs(rho_calc - rho0))
  if (emax > epsilon & errorloop == TRUE) {
    EL <- error_loop(k_cat = k_cat, k_cont = k_cont, k_pois = k_pois,
                     k_nb = k_nb, Y_cat = Y_cat, Y = Y, Yb = Yb,
                     Y_pois = Y_pois, Y_nb = Y_nb, marginal = marginal,
                     support = support, method = method, means = means,
                     vars = vars, constants = constants,
                     lam = lam, size = size, prob = prob, mu = mu,
                     n = n, seed = seed, epsilon = epsilon, maxit = maxit,
                     rho0 = rho0, Sigma = Sigma, rho_calc = rho_calc,
                     extra_correct = extra_correct)
    Sigma <- EL$Sigma
    rho_calc <- EL$rho_calc
    Y_cat <- EL$Y_cat
    Y <- EL$Y
    Yb <- EL$Yb
    Y_pois <- EL$Y_pois
    Y_nb <- EL$Y_nb
    niter <- EL$niter
  }
  stop.time.error <- Sys.time()
  Sigma2 <- Sigma
  emax <- max(abs(rho_calc - rho0))
  niter <- as.data.frame(niter)
  rownames(niter) <- c(1:k)
  colnames(niter) <- c(1:k)
  stop.time <- Sys.time()
  Time.constants <- round(difftime(stop.time.constants, start.time.constants,
                                   units = "min"), 3)
  cat("\nConstants calculation time:", Time.constants, "minutes \n")
  Time.intercorr <- round(difftime(stop.time.intercorr, start.time.intercorr,
                                   units = "min"), 3)
  cat("Intercorrelation calculation time:", Time.intercorr, "minutes \n")
  Time.error <- round(difftime(stop.time.error, start.time.error,
                               units = "min"), 3)
  cat("Error loop calculation time:", Time.error, "minutes \n")
  Time <- round(difftime(stop.time, start.time, units = "min"), 3)
  cat("Total Simulation time:", Time, "minutes \n")
  result <- list()
  if (k_cat > 0) {
    summary_cat <- list()
    for (i in 1:k_cat) {
      summary_cat[[i]] <- as.data.frame(cbind(append(marginal[[i]], 1),
                                              cumsum(table(Y_cat[, i]))/n))
      colnames(summary_cat[[i]]) = c("Target", "Simulated")
    }
    result <- append(result, list(ordinal_variables = as.data.frame(Y_cat),
                                  summary_ordinal = summary_cat))
  }
  if (k_cont > 0) {
    cont_sum <- describe(Yb, type = 1)
    sim_fifths <- rep(NA, k_cont)
    sim_sixths <- rep(NA, k_cont)
    for (i in 1:k_cont) {
      sim_fifths[i] <- calc_moments(Yb[, i])[5]
      sim_sixths[i] <- calc_moments(Yb[, i])[6]
    }
    cont_sum <- as.data.frame(cbind(c(1:k_cont),
                                    cont_sum[, -c(1, 6, 7, 10, 13)],
                                    sim_fifths, sim_sixths))
    colnames(cont_sum) <- c("Distribution", "n", "mean", "sd", "median",
                            "min", "max", "skew", "skurtosis", "fifth",
                            "sixth")
    target_sum <- as.data.frame(cbind(c(1:k_cont), means, sqrt(vars), skews,
                                      skurts, fifths, sixths))
    colnames(target_sum) <- c("Distribution", "mean", "sd", "skew",
                              "skurtosis", "fifth", "sixth")
    result <- append(result, list(constants = as.data.frame(constants),
                                  continuous_variables = as.data.frame(Yb),
                                  summary_continuous = cont_sum,
                                  summary_targetcont = target_sum,
                                  sixth_correction = SixCorr,
                                  valid.pdf = Valid.PDF))
  }
  if (k_pois > 0) {
    summary_pois <- describe(Y_pois, type = 1)
    summary_pois <- as.data.frame(cbind(summary_pois$vars, summary_pois$n,
                                      summary_pois$mean, lam,
                                      (summary_pois[, 4])^2, lam,
                                      summary_pois$median, summary_pois$min,
                                      summary_pois$max, summary_pois$skew,
                                      summary_pois$kurtosis))
    colnames(summary_pois) <- c("Distribution", "n", "mean", "Exp_mean",
                                "var", "Exp_var", "median", "min", "max",
                                "skew", "skurtosis")
    result <- append(result, list(Poisson_variables = as.data.frame(Y_pois),
                                  summary_Poisson = summary_pois))
  }
  if (k_nb > 0) {
    summary_nb <- describe(Y_nb, type = 1)
    if (length(prob) > 0) {
      mu <- size * (1 - prob)/prob
    }
    if (length(mu) > 0) {
      prob <- size/(mu + size)
    }
    summary_nb <- as.data.frame(cbind(summary_nb$vars, summary_nb$n, prob,
                                      summary_nb$mean, mu,
                                      (summary_nb[, 4])^2,
                                      size * (1 - prob)/prob^2,
                                      summary_nb$median, summary_nb$min,
                                      summary_nb$max,
                                      summary_nb$skew, summary_nb$kurtosis))
    colnames(summary_nb) <- c("Distribution", "n", "success_prob", "mean",
                              "Exp_mean", "var", "Exp_var", "median",
                              "min", "max", "skew", "skurtosis")
    result <- append(result, list(Neg_Bin_variables = as.data.frame(Y_nb),
                                  summary_Neg_Bin = summary_nb))
  }
  result <- append(result, list(correlations = rho_calc, Sigma1 = Sigma1,
                                Sigma2 = Sigma2,
                                Constants_Time = Time.constants,
                                Intercorrelation_Time = Time.intercorr,
                                Error_Loop_Time = Time.error,
                                Simulation_Time = Time,
                                niter = niter, maxerr = emax))
  result
}
