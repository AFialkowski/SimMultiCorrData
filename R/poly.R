#' @title Headrick's Fifth-Order Polynomial Transformation Equations
#'
#' @description This function contains Headrick's fifth-order polynomial transformation equations
#'     (2002, \doi{10.1016/S0167-9473(02)00072-5}).  It is used in
#'     \code{\link[SimMultiCorrData]{find_constants}} to find the constants c1, c2, c3, c4, and c5 (\eqn{c0 = -c2 - 3 * c4})
#'     that satisfy the equations given skewness, standardized kurtosis, and standardized fifth and sixth cumulant values.
#'     It can be used to verify a set of constants satisfy the equations.  Note that there exist solutions that yield
#'     invalid power method pdfs (see \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{pdf_check}}).  This function would not ordinarily be called by the user.
#' @param c a vector of constants c1, c2, c3, c4, c5; note that \code{\link[SimMultiCorrData]{find_constants}} returns
#'     c0, c1, c2, c3, c4, c5
#' @param a a vector c(skewness, standardized kurtosis, standardized fifth cumulant, standardized sixth cumulant)
#' @export
#' @keywords constants, Headrick
#' @seealso \code{\link[SimMultiCorrData]{fleish}}, \code{\link[SimMultiCorrData]{power_norm_corr}},
#'     \code{\link[SimMultiCorrData]{pdf_check}}, \code{\link[SimMultiCorrData]{find_constants}}
#' @return a list of length 5; if the constants satisfy the equations, returns 0 for all list elements
#' @references
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
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3), 1 - 17. \doi{10.18637/jss.v019.i03}.
#'
#' @examples
#' # Laplace Distribution
#' poly(c = c(0.727709, 0, 0.096303, 0, -0.002232), a = c(0, 3, 0, 30))
poly <- function(c, a) {
  f <- numeric(5)
  g1 <- a[1]
  g2 <- a[2]
  g3 <- a[3]
  g4 <- a[4]
  f[1] <- (c[1]^2 + 2 * c[2]^2 + 24 * c[2] * c[4] + 6 * c[1] *
            (c[3] + 5 * c[5]) + 3 * (5 * c[3]^2 + 32 * c[4]^2 +
                                       70 * c[3] * c[5] +
                                       315 * c[5]^2)) - 1
  f[2] <- (2 * (4 * c[2]^3 + 108 * c[2]^2 * c[4] + 3 * c[1]^2 *
                 (c[2] + 6 * c[4]) + 18 * c[1] * (2 * c[2] * c[3] +
                                                    16 * c[3] * c[4] +
                                                    15 * c[2] * c[5] +
                                                    150 * c[4] * c[5]) +
                 9 * c[2] * (15 * c[3]^2 + 128 * c[4]^2 +
                               280 * c[3] * c[5] + 1575 * c[5]^2) +
                 54 * c[4] * (25 * c[3]^2 + 88 * c[4]^2 +
                                560 * c[3] * c[5] + 3675 * c[5]^2))) - g1
  f[3] <- (24 * (2 * c[2]^4 + 96 * c[2]^3 * c[4] + c[1]^3 * (c[3] + 10 *
                                                               c[5]) +
                  30 * c[2]^2 * (6 * c[3]^2 + 64 * c[4]^2 + 140 * c[3] * c[5] +
                                   945 * c[5]^2) +
                  c[1]^2 * (2 * c[2]^2 + 18 * c[3]^2 + 36 * c[2] * c[4] +
                              192 * c[4]^2 + 375 * c[3] * c[5] +
                              2250 * c[5]^2) +
                  36 * c[2] * c[4] * (125 * c[3]^2 + 528 * c[4]^2 +
                                        3360 * c[3] * c[5] + 25725 * c[5]^2) +
                  3 * c[1] * (45 * c[3]^3 + 1584 * c[3] * c[4]^2 +
                                1590 * c[3]^2 * c[5] + 21360 * c[4]^2 * c[5] +
                                21525 * c[3] * c[5]^2 + 110250 * c[5]^3 +
                                12 * c[2]^2 * (c[3] + 10 * c[5]) +
                                8 * c[2] * c[4] * (32 * c[3] + 375 * c[5])) +
                  9 * (45 * c[3]^4 + 8704 * c[4]^4 + 2415 * c[3]^3 * c[5] +
                         932400 * c[4]^2 * c[5]^2 + 3018750 * c[5]^4 +
                         20 * c[3]^2 * (178 * c[4]^2 + 2765 * c[5]^2) +
                         35 * c[3] * (3104 * c[4]^2 * c[5] +
                                        18075 * c[5]^3)))) - g2
  f[4] <- (24 * (16 * c[2]^5 + 5 * c[1]^4 * c[4] + 1200 * c[2]^4 * c[4] +
                  10 * c[1]^3 * (3 * c[2] * c[3] + 42 * c[3] * c[4] +
                                   40 * c[2] * c[5] + 570 * c[4] * c[5]) +
                  300 * c[2]^3 * (10 * c[3]^2 + 128 * c[4]^2 +
                                    280 * c[3] * c[5] + 2205 * c[5]^2) +
                  1080 * c[2]^2 * c[4] *
                  (125 * c[3]^2 + 3920 * c[3] * c[5] +
                     28 * (22 * c[4]^2 + 1225 * c[5]^2)) +
                  10 * c[1]^2 * (2 * c[2]^3 + 72 * c[2]^2 * c[4] +
                                   3 * c[2] * (24 * c[3]^2 + 320 * c[4]^2 +
                                                 625 * c[3] * c[5] +
                                                 4500 * c[5]^2) +
                                   9 * c[4] * (109 * c[3]^2 +
                                                 528 * c[4]^2 +
                                                 3130 * c[3] * c[5] +
                                                 24975 * c[5]^2)) +
                  30 * c[1] * (8 * c[2]^3 * (2 * c[3] + 25 * c[5]) +
                                 40 * c[2]^2 * c[4] * (16 * c[3] +
                                                         225 * c[5]) +
                                 3 * c[2] * (75 * c[3]^3 +
                                               3168 * c[3] * c[4]^2 +
                                               3180 * c[3]^2 * c[5] +
                                               49840 * c[4]^2 * c[5] +
                                               50225 * c[3] * c[5]^2 +
                                               294000 * c[5]^3) +
                                 6 * c[4] * (555 * c[3]^3 +
                                               8704 * c[3] * c[4]^2 +
                                               26225 * c[3]^2 * c[5] +
                                               152160 * c[4]^2 * c[5] +
                                               459375 * c[3] * c[5]^2 +
                                               2963625 * c[5]^3)) +
                  90 * c[2] * (270 * c[3]^4 + 16905 * c[3]^3 * c[5] +
                                 280 * c[3]^2 * (89 * c[4]^2 +
                                                   1580 * c[5]^2) +
                                 35 * c[3] * (24832 * c[4]^2 * c[5] +
                                                162675 * c[5]^3) +
                                 4 * (17408 * c[4]^4 +
                                        2097900 * c[4]^2 * c[5]^2 +
                                        7546875 * c[5]^4)) +
                  27 * c[4] * (14775 * c[3]^4 + 1028300 * c[3]^3 * c[5] +
                                 50 * c[3]^2 * (10144 * c[4]^2 +
                                                  594055 * c[5]^2) +
                                 700 * c[3] * (27904 * c[4]^2 * c[5] +
                                                 598575 * c[5]^3) +
                                 3 * (316928 * c[4]^4 +
                                        68908000 * c[4]^2 * c[5]^2 +
                                        806378125 * c[5]^4)))) - g3
  f[5] <- (120 * (32 * c[2]^6 + 3456 * c[2]^5 * c[4] + 6 * c[1]^5 * c[5] +
                   3 * c[1]^4 * (9 * c[3]^2 + 16 * c[2] * c[4] + 168 * c[4]^2 +
                                   330 * c[3] * c[5] + 2850 * c[5]^2) +
                   720 * c[2]^4 * (15 * c[3]^2 + 224 * c[4]^2 +
                                     490 * c[3] * c[5] + 4410 * c[5]^2) +
                   6048 * c[2]^3 * c[4] * (125 * c[3]^2 + 704 * c[4]^2 +
                                             4480 * c[3] * c[5] +
                                             44100 * c[5]^2) +
                   12 * c[1]^3 * (4 * c[2]^2 * (3 * c[3] + 50 * c[5]) +
                                    60 * c[2] * c[4] * (7 * c[3] + 114 *
                                                          c[5]) +
                                    3 * (24 * c[3]^3 + 1192 * c[3] * c[4]^2 +
                                           1170 * c[3]^2 * c[5] +
                                           20440 * c[4]^2 * c[5] +
                                           20150 * c[3] * c[5]^2 +
                                           124875 * c[5]^3)) +
                   216 * c[2]^2 * (945 * c[3]^4 + 67620 * c[3]^3 * c[5] +
                                     560 * c[3]^2 * (178 * c[4]^2 +
                                                       3555 * c[5]^2) +
                                     315 * c[3] * (12416 * c[4]^2 * c[5] +
                                                     90375 * c[5]^3) +
                                     6 * (52224 * c[4]^4 +
                                            6993000 * c[4]^2 * c[5]^2 +
                                            27671875 * c[5]^4)) +
                   6 * c[1]^2 * (8 * c[2]^4 + 480 * c[2]^3 * c[4] +
                                   180 * c[2]^2 * (4 * c[3]^2 + 64 * c[4]^2 +
                                                     125 * c[3] * c[5] +
                                                     1050 * c[5]^2) +
                                   72 * c[2] * c[4] * (327 * c[3]^2 +
                                                         1848 * c[4]^2 +
                                                         10955 * c[3] * c[5] +
                                                         99900 * c[5]^2) +
                                   9 * (225 * c[3]^4 + 22824 * c[3]^2 * c[4]^2 +
                                          69632 * c[4]^4 + 15090 * c[3]^3 *
                                          c[5] +
                                          830240 * c[3] * c[4]^2 * c[5] +
                                          412925 * c[3]^2 * c[5]^2 +
                                          8239800 * c[4]^2 * c[5]^2 +
                                          5475750 * c[3] * c[5]^3 +
                                          29636250 * c[5]^4)) +
                   1296 * c[2] * c[4] * (5910 * c[3]^4 + 462735 * c[3]^3 *
                                           c[5] +
                                           c[3]^2 * (228240 * c[4]^2 +
                                                       14851375 * c[5]^2) +
                                           175 * c[3] * (55808 * c[4]^2 *
                                                           c[5] +
                                                           1316865 * c[5]^3) +
                                           3 * (158464 * c[4]^4 +
                                                  37899400 * c[4]^2 * c[5]^2 +
                                                  483826875 * c[5]^4)) +
                   27 * (9945 * c[3]^6 + 92930048 * c[4]^6 +
                           1166130 * c[3]^5 * c[5] + 35724729600 * c[4]^4 *
                           c[5]^2 +
                           977816385000 * c[4]^2 * c[5]^4 + 1907724656250 *
                           c[5]^6 +
                           180 * c[3]^4 * (16082 * c[4]^2 + 345905 * c[5]^2) +
                           140 * c[3]^3 * (1765608 * c[4]^2 * c[5] +
                                             13775375 * c[5]^3) +
                           15 * c[3]^2 * (4076032 * c[4]^4 +
                                            574146160 * c[4]^2 * c[5]^2 +
                                            2424667875 * c[5]^4) +
                           210 * c[3] * (13526272 * c[4]^4 * c[5] +
                                           687499200 * c[4]^2 * c[5]^3 +
                                           1876468125 * c[5]^5)) +
                   18 * c[1] * (80 * c[2]^4 * (c[3] + 15 * c[5]) +
                                  160 * c[2]^3 * c[4] * (32 * c[3] + 525 *
                                                           c[5]) +
                                  12 * c[2]^2 *
                                  (225 * c[3]^3 + 11088 * c[3] * c[4]^2 +
                                     11130 * c[3]^2 * c[5] +
                                     199360 * c[4]^2 * c[5] +
                                     200900 * c[3] * c[5]^2 +
                                     1323000 * c[5]^3) +
                                  24 * c[2] * c[4] *
                                  (3885 * c[3]^3 + 69632 * c[3] * c[4]^2 +
                                     209800 * c[3]^2 * c[5] +
                                     1369440 * c[4]^2 * c[5] +
                                     4134375 * c[3] * c[5]^2 + 29636250 *
                                     c[5]^3) +
                                  9 * (540 * c[3]^5 + 48585 * c[3]^4 * c[5] +
                                         20 * c[3]^3 * (4856 * c[4]^2 +
                                                          95655 * c[5]^2) +
                                         80 * c[3]^2 * (71597 * c[4]^2 * c[5] +
                                                          513625 * c[5]^3) +
                                         4 * c[3] * (237696 * c[4]^4 +
                                                       30726500 * c[4]^2 *
                                                       c[5]^2 +
                                                       119844375 * c[5]^4) +
                                         5 * c[5] * (4076032 * c[4]^4 +
                                                       191074800 * c[4]^2 *
                                                       c[5]^2 +
                                                       483826875 *
                                                       c[5]^4))))) - g4
  return(f)
}
