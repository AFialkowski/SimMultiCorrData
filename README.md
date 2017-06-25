# SimMultiCorrData
The goal of `SimMultiCorrData` is to generate continuous (normal or non-normal), binary, ordinal, and count (Poisson or Negative Binomial) variables with a specified correlation matrix.  It can also produce a single continuous variable.  This package can be used to simulate data sets that mimic real-world situations (i.e. clinical data sets, plasmodes, as in Vaughan et al., 2009).  All variables are generated from standard normal variables with an imposed intermediate correlation matrix.  Continuous variables are simulated by specifying mean, variance, skewness, standardized kurtosis, and fifth and sixth standardized cumulants using either Fleishman's Third-Order or Headrick's Fifth-Order Polynomial Transformation.  Binary and ordinal variables are simulated using a modification of `GenOrd::ordsample`.  Count variables are simulated using the inverse cdf method.  There are two simulation pathways which differ primarily according to the calculation of the intermediate correlation matrix `Sigma`.  In **Method 1**, the intercorrelations involving count variables are determined using a simulation based, logarithmic correlation correction (adapting Yahav and Shmueli's 2012 method).  In **Method 2**, the count variables are treated as ordinal (adapting Barbiero and Ferrari's 2015 modification of `GenOrd`).  There is an optional error loop that corrects the final correlation matrix to be within a user-specified precision value. The package also includes functions to calculate standardized cumulants for theoretical distributions or from real data sets, check if a target correlation matrix is within the possible correlation bounds (given the distributions of the simulated variables), summarize results (numerically or graphically), to verify valid power method pdfs, and to calculate lower standardized kurtosis bounds.

There are several vignettes which accompany this package that may help the user understand the simulation and analysis methods.

1) **Benefits of SimMultiCorrData and Comparison to Other Packages** describes some of the ways `SimMultiCorrData` improves
upon other simulation packages.

2) **Variable Types** describes the different types of variables that can be simulated in `SimMultiCorrData`.

3) **Function by Topic** describes each function, separated by topic.

4) **Comparison of Method 1 and Method 2** describes the two simulation pathways that can be followed.

5) **Overview of Error Loop** details the algorithm involved in the optional error loop that improves the accuracy of the
simulated variables' correlation matrix.

6) **Overall Workflow for Data Simulation** gives a step-by-step guideline to follow with an example containing continuous
(normal and non-normal), binary, ordinal, Poisson, and Negative Binomial variables.  It also demonstrates the use of the
standardized cumulant calculation function, correlation check functions, the lower kurtosis boundary function, and the plotting functions.

7) **Comparison of Simulation Distribution to Theoretical Distribution or Empirical Data** gives a step-by-step guideline for
comparing a simulated univariate continuous distribution to the target distribution with an example.

8) **Using the Sixth Cumulant Correction to Find Valid Power Method Pdfs** demonstrates how to use the sixth cumulant correction
to generate a valid power method pdf and the effects this has on the resulting distribution.
