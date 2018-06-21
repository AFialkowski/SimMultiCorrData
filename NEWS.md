# SimMultiCorrData 0.2.2      
1. Fixed error in formation of intermediate matrices with `findintercorr2`.     
2. Fixed check for repeated continuous distribution when using `method` = "Fleishman".
3. Updated citations.
4. Fixed function/package links in documentation.
5. Added automated testing with 'testthat' package.

# SimMultiCorrData 0.2.1
1. Updated `calc_theory()` and plotting functions which call it to permit *pdf* specified by `fx`, `lower`, and `upper`.
2. Adjusted `rcorrvar()` and `rcorrvar2()` summary of continuous variables when using `method = "Fleishman"`.
3. Fixed error in `rcorrvar()`, `rcorrvar2()`, `valid_corr()`, `valid_corr2()`, and `error_loop()` to permit 0 or 1 continuous variables.
4. Updated `calc_lower_skurt()` for case of non-convergence when applying `Six` vector with `method = "Polynomial"`.
5. Added `lower` and `upper` parameters to `plot_cdf()` to use as inputs for `cdf_prob()`.

# SimMultiCorrData 0.2.0
1. Fixed error in `findintercorr2()` so now you can generate 1 ordinal variable using correlation method 2 (with `rcorrvar2()`).
2. Fixed error in `chat_nb()` so you can use `size` (success probability) and `mu` (mean) parameters for Negative Binomial variables when using correlation method 1 (with `rcorrvar1()`).
3. Updated `find_constants()` and `calc_lower_skurt()` (to remove duplicate rows in solutions before executing `pdf_check()`) in order to decrease computation time.
4. Updated `rcorrvar()`, `rcorrvar2()`, `valid_corr()`, and `valid_corr2()` to check for identical continuous distributions before calculating the power method constants in order to decrease computation time.  If a distribution is repeated, the constants are only calculated once.
5. Made the following updates to `error_loop()` and `error_vars()`:
  * all variables are regenerated in each iteration and the final correlation is calculated each time and at the end with all variables (instead of pairwise)
  * eigenvalue decomposition on `Sigma` is done using the maximum of 0 and the eigenvalues (in case `Sigma` is not positive-definite and the eigenvalues are negative); this replaces the use of `Matrix::nearPD()`
  * update function is based on the correlation calculated from the previous iteration instead of the pre-error loop final correlation.
  * fixed `ifelse()` statement in choice of update function (affects negative correlations only)
6. Updated the **Overview of Error Loop** vignette to reflect above changes.
7. Fixed `ifelse()` statement in choice of update function for `ordnorm()` (affects negative correlations only).
8. Made the following updates to `calc_theory()`:
 * `params` input accepts up to 4 parameters
 *  39 distributions are available by name (as `Dist` input)
9. Made the following updates to `plot_pdf_theory()`, `plot_sim_pdf_theory()`, and `plot_sim_theory()`:
 * `params` input accepts up to 4 parameters
 *  39 distributions are available by name (as `Dist` input) plus Poisson and Negative Binomial for `plot_sim_pdf_theory()` and `plot_sim_theory()`
10. Added extra `ggplot2` parameters to the graphing functions to allow control over the appearance of the legend, axes labels and titles, and plot title.
11. Changed the example in the **Overall Workflow for Data Simulation** vignette.
12. Changed the examples in the `rcorrvar()` and `rcorrvar2()` documentation.
13. Updated documentation to some of the functions.


# SimMultiCorrData 0.1.0
Initial package release.
