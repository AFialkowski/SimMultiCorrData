#' @title Calculate Intermediate MVN Correlation for Ordinal, Continuous, Poisson, or Negative Binomial Variables: Method 1
#'
#' @description This function calculates a \code{k x k} intermediate matrix of correlations, where \code{k = k_cat + k_cont + k_pois + k_nb},
#'     to be used in simulating variables with \code{\link[SimMultiCorrData]{rcorrvar}}.  The ordering of the variables must be
#'     ordinal, continuous, Poisson, and Negative Binomial (note that it is possible for \code{k_cat}, \code{k_cont}, \code{k_pois},
#'     and/or \code{k_nb} to be 0).
#'     The function first checks that the target correlation matrix \code{rho} is positive-definite and the marginal distributions for the
#'     ordinal variables are cumulative probabilities with r - 1 values (for r categories).  There is a warning given at the end of simulation
#'     if the calculated intermediate correlation matrix \code{Sigma} is not positive-definite.  This function is called by the simulation function
#'     \code{\link[SimMultiCorrData]{rcorrvar}}, and would only be used separately if the user wants to find the intermediate correlation matrix
#'     only.  The simulation functions also return the intermediate correlation matrix.
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
#'     The processes used to find the intermediate correlations for each variable type are described below.  Please see the
#'     corresponding function help page for more information:
#'
#' @section Ordinal Variables:
#' Correlations are computed pairwise.  If both variables are binary, the method of Demirtas et al. (2012) is used to find the
#' tetrachoric correlation (code adapted from \code{\link[BinNonNor]{Tetra.Corr.BB}}).  This method is based on Emrich and Piedmonte's
#' (1991) work, in which the joint binary distribution is determined from the third and higher moments of a multivariate normal
#' distribution: Let \eqn{Y_{1}} and \eqn{Y_{2}} be binary variables with \eqn{E[Y_{1}] = Pr(Y_{1} = 1) = p_{1}},
#' \eqn{E[Y_{2}] = Pr(Y_{2} = 1) = p_{2}}, and correlation \eqn{\rho_{y1y2}}.  Let \eqn{\Phi[x_{1}, x_{2}, \rho_{x1x2}]} be the
#' standard bivariate normal cumulative distribution function, given by:
#' \deqn{\Phi[x_{1}, x_{2}, \rho_{x1x2}] = \int_{-\infty}^{x_{1}} \int_{-\infty}^{x_{2}} f(z_{1}, z_{2}, \rho_{x1x2}) dz_{1} dz_{2}}
#' where
#' \deqn{f(z_{1}, z_{2}, \rho_{x1x2}) = [2\pi\sqrt{1 - \rho_{x1x2}^2}]^{-1} * exp[-0.5(z_{1}^2 - 2\rho_{x1x2}z_{1}z_{2} + z_{2}^2)/(1 - \rho_{x1x2}^2)]}
#' Then solving the equation
#' \deqn{\Phi[z(p_{1}), z(p_{2}), \rho_{x1x2}] = \rho_{y1y2}\sqrt{p_{1}(1 - p_{1})p_{2}(1 - p_{2})} + p_{1}p_{2}}
#' for \eqn{\rho_{x1x2}} gives the intermediate correlation of the standard normal variables needed to generate binary variables with
#' correlation \eqn{\rho_{y1y2}}.  Here \eqn{z(p)} indicates the \eqn{pth} quantile of the standard normal distribution.
#'
#' Otherwise, \code{\link[SimMultiCorrData]{ordnorm}} is called for each pair.  If the resulting
#' intermediate matrix is not positive-definite, there is a warning given because it may not be possible to find a MVN correlation
#' matrix that will produce the desired marginal distributions after discretization.  Any problems with positive-definiteness are
#' corrected later.
#'
#' @section Continuous Variables:
#' Correlations are computed pairwise.  \code{\link[SimMultiCorrData]{findintercorr_cont}} is called for each pair.
#'
#' @section Poisson Variables:
#' \code{\link[SimMultiCorrData]{findintercorr_pois}} is called to calculate the intermediate MVN correlation for all Poisson variables.
#'
#' @section Negative Binomial Variables:
#' \code{\link[SimMultiCorrData]{findintercorr_nb}} is called to calculate the intermediate MVN correlation for all Negative
#' Binomial variables.
#'
#' @section Continuous - Ordinal Pairs:
#' \code{\link[SimMultiCorrData]{findintercorr_cont_cat}} is called to calculate the intermediate MVN correlation for all
#' Continuous and Ordinal combinations.
#'
#' @section Ordinal - Poisson Pairs:
#' \code{\link[SimMultiCorrData]{findintercorr_cat_pois}} is called to calculate the intermediate MVN correlation for all
#' Ordinal and Poisson combinations.
#'
#' @section Ordinal - Negative Binomial Pairs:
#' \code{\link[SimMultiCorrData]{findintercorr_cat_nb}} is called to calculate the intermediate MVN correlation for all
#' Ordinal and Negative Binomial combinations.
#'
#' @section Continuous - Poisson Pairs:
#' \code{\link[SimMultiCorrData]{findintercorr_cont_pois}} is called to calculate the intermediate MVN correlation for all
#' Continuous and Poisson combinations.
#'
#' @section Continuous - Negative Binomial Pairs:
#' \code{\link[SimMultiCorrData]{findintercorr_cont_nb}} is called to calculate the intermediate MVN correlation for all
#' Continuous and Negative Binomial combinations.
#'
#' @section Poisson - Negative Binomial Pairs:
#' \code{\link[SimMultiCorrData]{findintercorr_pois_nb}} is called to calculate the intermediate MVN correlation for all
#' Poisson and Negative Binomial combinations.
#'
#' @param n the sample size (i.e. the length of each simulated variable)
#' @param k_cont the number of continuous variables (default = 0)
#' @param k_cat the number of ordinal (r >= 2 categories) variables (default = 0)
#' @param k_pois the number of Poisson variables (default = 0)
#' @param k_nb the number of Negative Binomial variables (default = 0)
#' @param method the method used to generate the \code{k_cont} continuous variables.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param constants a matrix with \code{k_cont} rows, each a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or
#'     c0, c1, c2, c3, c4, c5 (if \code{method} = "Polynomial") like that returned by
#'     \code{\link[SimMultiCorrData]{find_constants}}
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1; default = list())
#' @param support a list of length equal to \code{k_cat}; the i-th element is a vector of containing the r
#'     ordered support values; if not provided (i.e. \code{support} = list()), the default is for the i-th element to be the vector 1, ..., r
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{dpois}})
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}})
#' @param prob a vector of success probability parameters
#' @param mu a vector of mean parameters (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL)
#' @param rho the target correlation matrix (\emph{must be ordered ordinal, continuous, Poisson, Negative Binomial}; default = NULL)
#' @param seed the seed value for random number generation (default = 1234)
#' @param epsilon the maximum acceptable error between the final and target correlation matrices (default = 0.001)
#'     in the calculation of ordinal intermediate correlations with \code{\link[SimMultiCorrData]{ordnorm}}
#' @param maxit the maximum number of iterations to use (default = 1000) in the calculation of ordinal
#'     intermediate correlations with \code{\link[SimMultiCorrData]{ordnorm}}
#' @importFrom psych describe
#' @import stats
#' @import utils
#' @import BB
#' @importFrom Matrix nearPD
#' @import GenOrd
#' @export
#' @keywords intermediate, correlation, continuous, ordinal, Poisson, Negative Binomial, Fleishman, Headrick, method1
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimMultiCorrData]{rcorrvar}}
#' @return the intermediate MVN correlation matrix
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
#' Inan G & Demirtas H (2016). BinNonNor: Data Generation with Binary and Continuous Non-Normal Components.
#'     R package version 1.3. \url{https://CRAN.R-project.org/package=BinNonNor}
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
#' @examples \dontrun{
#'
#' # Binary, Ordinal, Continuous, Poisson, and Negative Binomial Variables
#'
#' options(scipen = 999)
#' seed <- 1234
#' n <- 10000
#'
#' # Continuous Distributions: Normal, t (df = 10), Chisq (df = 4),
#' # Beta (a = 4, b = 2), Gamma (a = 4, b = 4)
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
#' # calculate constants
#' con <- matrix(1, nrow = ncol(M), ncol = 6)
#' for (i in 1:ncol(M)) {
#'  con[i, ] <- find_constants(method = "Polynomial", skews = M[1, i],
#'                             skurts = M[2, i], fifths = M[3, i],
#'                             sixths = M[4, i], Six = NULL,
#'                             cstart = NULL, n = 25, seed = seed)
#' }
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
#'                     fifths = M[3, ], sixths = M[4, ], Six = NULL,
#'                     marginal = marginal, lam = lam, size = size,
#'                     prob = prob, mu = NULL, rho = Rey, n = 100000,
#'                     seed = seed)
#'
#' # Find intermediate correlation
#' Sigma1 <- findintercorr(n = n, k_cont = ncont, k_cat = ncat, k_pois = npois,
#'                         k_nb = nnb, method = "Polynomial", constants = con,
#'                         marginal = marginal, support = list(),
#'                         nrand = 100000, lam = lam, size = size, prob = prob,
#'                         mu = NULL, rho = Rey, seed = seed, epsilon = 0.001,
#'                         maxit = 1000)
#' Sigma1
#'
#' }
findintercorr <- function(n, k_cont = 0, k_cat = 0, k_pois = 0, k_nb = 0,
                          method = c("Fleishman", "Polynomial"), constants,
                          marginal = list(), support = list(), nrand = 100000,
                          lam = NULL, size = NULL, prob = NULL, mu = NULL,
                          rho = NULL, seed = 1234, epsilon = 0.001,
                          maxit = 1000) {
  k <- k_cat + k_cont + k_pois + k_nb
  Spearman <- FALSE
  if (ncol(rho) != k) {
    stop("Dimension of correlation matrix does not match the number of
         variables!")
  }
  if (!isSymmetric(rho) | min(eigen(rho, symmetric = TRUE)$values) < 0 |
      !all(diag(rho) == 1)) {
    stop("Correlation matrix not valid!")
  }
  if (k_cat > 0) {
    if (!all(unlist(lapply(marginal,
                           function(x) (sort(x) == x & min(x) > 0 &
                                        max(x) < 1))))) {
      stop("Error with given marginal distributions!")
    }
    len <- length(support)
    kj <- numeric(k_cat)
    for (i in 1:k_cat) {
      kj[i] <- length(marginal[[i]]) + 1
      if (len == 0) {
        support[[i]] <- 1:kj[i]
      }
    }
  }
  rho_list <- separate_rho(k_cat = k_cat, k_cont = k_cont, k_pois = k_pois,
                           k_nb = k_nb, rho = rho)
  rho_cat <- rho_list$rho_cat
  rho_cat_pois <- rho_list$rho_cat_pois
  rho_cat_nb <- rho_list$rho_cat_nb
  rho_cont_cat <- rho_list$rho_cont_cat
  rho_cont <- rho_list$rho_cont
  rho_cont_pois <- rho_list$rho_cont_pois
  rho_cont_nb <- rho_list$rho_cont_nb
  rho_pois <- rho_list$rho_pois
  rho_pois_nb <- rho_list$rho_pois_nb
  rho_nb <- rho_list$rho_nb
  Sigma_cat <- NA
  Sigma_cat_cont <- NA
  Sigma_cat_pois <- NA
  Sigma_cat_nb <- NA
  Sigma_cont_cat <- NA
  Sigma_cont <- NA
  Sigma_cont_pois <- NA
  Sigma_cont_nb <- NA
  Sigma_pois_cat <- NA
  Sigma_pois_cont <- NA
  Sigma_pois <- NA
  Sigma_pois_nb <- NA
  Sigma_nb_cat <- NA
  Sigma_nb_cont <- NA
  Sigma_nb_pois <- NA
  Sigma_nb <- NA
  if (k_cat == 1) {
    Sigma_cat <- matrix(1, nrow = k_cat, ncol = k_cat)
  }
  if (k_cat > 1) {
    Sigma_cat <- diag(1, k_cat, k_cat)
    for (i in 1:(k_cat - 1)) {
      for (j in (i + 1):k_cat) {
        if (length(marginal[[i]]) == 1 & length(marginal[[j]]) == 1) {
          corr_bin <- function(rho) {
            phix1x2 <- integrate(function(z2) {
              sapply(z2, function(z2) {
                integrate(function(z1) ((2 * pi * sqrt((1 - rho^2)))^-1) *
                    exp(-(z1^2 - 2 * rho * z1 * z2 + z2^2)/(2 * (1 - rho^2))),
                    -Inf, qnorm(1 - marginal[[i]][1]))$value
              })
            }, -Inf, qnorm(1 - marginal[[j]][1]))$value -
              rho_cat[i, j] * sqrt(marginal[[i]][1] * (1 - marginal[[j]][1]) *
                                  marginal[[j]][1] * (1 - marginal[[i]][1])) -
              ((1 - marginal[[i]][1]) * (1 - marginal[[j]][1]))
            phix1x2
          }
          Sigma_cat[i, j] <- suppressWarnings(dfsane(par = 0, fn = corr_bin,
                                           control = list(trace = FALSE)))$par
         } else {
          Sigma_cat[i, j] <-
            suppressWarnings(ordnorm(list(marginal[[i]], marginal[[j]]),
                                     matrix(c(1, rho_cat[i, j],
                                              rho_cat[i, j], 1), 2, 2),
                                     list(support[[i]], support[[j]]),
                                     epsilon = epsilon,
                                     maxit = maxit)$SigmaC[1, 2])
        }
        Sigma_cat[j, i] <- Sigma_cat[i, j]
      }
    }
    if (min(eigen(Sigma_cat, symmetric = TRUE)$values) < 0) {
      warning("It is not possible to find a correlation matrix for MVN
              ensuring rho for the ordinal variables.
              Try the error loop.")
    }
  }
  if (k_cont > 0) {
    Sigma_cont <- diag(1, k_cont, k_cont)
    for (i in 1:k_cont) {
      for (j in (1:k_cont)) {
        if (j > i) {
          Sigma_cont[i, j] <-
            findintercorr_cont(method, constants[c(i, j), ],
                               matrix(rho_cont[i, j],
                                      nrow = 1, ncol = 1))$x
          Sigma_cont[j, i] <- Sigma_cont[i, j]
        }
      }
    }
  }
  if (k_cat > 0 & k_cont > 0) {
    Sigma_cont_cat <- findintercorr_cont_cat(method, constants,
                                             rho_cont_cat, marginal, support)
    Sigma_cat_cont <- t(Sigma_cont_cat)
  }
  if (k_cat > 0 & k_pois > 0) {
    Sigma_cat_pois <-
      findintercorr_cat_pois(rho_cat_pois, marginal = marginal,
                             lam, nrand, seed)
    Sigma_pois_cat <- t(Sigma_cat_pois)
  }
  if (k_cat > 0 & k_nb > 0) {
    Sigma_cat_nb <-
      findintercorr_cat_nb(rho_cat_nb, marginal = marginal,
                           size, prob, mu, nrand, seed)
    Sigma_nb_cat <- t(Sigma_cat_nb)
  }
  if (k_cont > 0 & k_pois > 0) {
    Sigma_cont_pois <- findintercorr_cont_pois(method, constants,
                                               rho_cont_pois, lam, nrand,
                                               seed)
    Sigma_pois_cont <- t(Sigma_cont_pois)
  }
  if (k_cont > 0 & k_nb > 0) {
    Sigma_cont_nb <-
      findintercorr_cont_nb(method, constants, rho_cont_nb, size,
                            prob, mu, nrand, seed)
    Sigma_nb_cont <- t(Sigma_cont_nb)
  }
  if (k_pois == 1) {
    Sigma_pois <- matrix(1, nrow = k_pois, ncol = k_pois)
  }
  if (k_pois > 1) {
    Sigma_pois <- findintercorr_pois(rho_pois, lam, nrand, seed)
  }
  if (k_nb == 1) {
    Sigma_nb <- matrix(1, nrow = k_nb, ncol = k_nb)
  }
  if (k_nb > 1) {
    Sigma_nb <- findintercorr_nb(rho_nb, size, prob, mu, nrand, seed)
  }
  if (k_pois > 0 & k_nb > 0) {
    Sigma_pois_nb <- findintercorr_pois_nb(rho_pois_nb, lam, size,
                                           prob, mu, nrand, seed)
    Sigma_nb_pois <- t(Sigma_pois_nb)
  }
  if (k_cat > 0 & k_cont == 0 & k_pois == 0 & k_nb == 0) {
    Sigma <- Sigma_cat
  }
  if (k_cat == 0 & k_cont > 0 & k_pois == 0 & k_nb == 0) {
    Sigma <- Sigma_cont
  }
  if (k_cat == 0 & k_cont == 0 & k_pois > 0 & k_nb == 0) {
    Sigma <- Sigma_pois
  }
  if (k_cat == 0 & k_cont == 0 & k_pois == 0 & k_nb > 0) {
    Sigma <- Sigma_nb
  }
  if (k_cat > 0 & k_cont > 0 & k_pois > 0 & k_nb > 0) {
    Sigma <- rbind(cbind(Sigma_cat, Sigma_cat_cont, Sigma_cat_pois,
                         Sigma_cat_nb),
                   cbind(Sigma_cont_cat, Sigma_cont, Sigma_cont_pois,
                         Sigma_cont_nb),
                   cbind(Sigma_pois_cat, Sigma_pois_cont, Sigma_pois,
                         Sigma_pois_nb),
                   cbind(Sigma_nb_cat, Sigma_nb_cont, Sigma_nb_pois,
                         Sigma_nb))
  }
  if (k_cat > 0 & k_cont > 0 & k_pois > 0 & k_nb == 0) {
    Sigma <- rbind(cbind(Sigma_cat, Sigma_cat_cont, Sigma_cat_pois),
                   cbind(Sigma_cont_cat, Sigma_cont, Sigma_cont_pois),
                   cbind(Sigma_pois_cat, Sigma_pois_cont, Sigma_pois))
  }
  if (k_cat > 0 & k_cont > 0 & k_pois == 0 & k_nb > 0) {
    Sigma <- rbind(cbind(Sigma_cat, Sigma_cat_cont, Sigma_cat_nb),
                   cbind(Sigma_cont_cat, Sigma_cont, Sigma_cont_nb),
                   cbind(Sigma_nb_cat, Sigma_nb_cont, Sigma_nb))
  }
  if (k_cat > 0 & k_cont == 0 & k_pois > 0 & k_nb > 0) {
    Sigma <- rbind(cbind(Sigma_cat, Sigma_cat_pois, Sigma_cat_nb),
                   cbind(Sigma_pois_cat, Sigma_pois, Sigma_pois_nb),
                   cbind(Sigma_nb_cat, Sigma_nb_pois, Sigma_nb))
  }
  if (k_cat == 0 & k_cont > 0 & k_pois > 0 & k_nb > 0) {
    Sigma <- rbind(cbind(Sigma_cont, Sigma_cont_pois, Sigma_cont_nb),
                   cbind(Sigma_pois_cont, Sigma_pois, Sigma_pois_nb),
                   cbind(Sigma_nb_cont, Sigma_nb_pois, Sigma_nb))
  }
  if (k_cat > 0 & k_cont > 0 & k_pois == 0 & k_nb == 0) {
    Sigma <- rbind(cbind(Sigma_cat, Sigma_cat_cont),
                   cbind(Sigma_cont_cat, Sigma_cont))
  }
  if (k_cat > 0 & k_cont == 0 & k_pois > 0 & k_nb == 0) {
    Sigma <- rbind(cbind(Sigma_cat, Sigma_cat_pois),
                   cbind(Sigma_pois_cat, Sigma_pois))
  }
  if (k_cat > 0 & k_cont == 0 & k_pois == 0 & k_nb > 0) {
    Sigma <- rbind(cbind(Sigma_cat, Sigma_cat_nb),
                   cbind(Sigma_nb_cat, Sigma_nb))
  }
  if (k_cat == 0 & k_cont > 0 & k_pois > 0 & k_nb == 0) {
    Sigma <- rbind(cbind(Sigma_cont, Sigma_cont_pois),
                   cbind(Sigma_pois_cont, Sigma_pois))
  }
  if (k_cat == 0 & k_cont > 0 & k_pois == 0 & k_nb > 0) {
    Sigma <- rbind(cbind(Sigma_cont, Sigma_cont_nb),
                   cbind(Sigma_nb_cont, Sigma_nb))
  }
  if (k_cat == 0 & k_cont == 0 & k_pois > 0 & k_nb > 0) {
    Sigma <- rbind(cbind(Sigma_pois, Sigma_pois_nb),
                   cbind(Sigma_nb_pois, Sigma_nb))
  }
  return(Sigma)
}
