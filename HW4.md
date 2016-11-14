# HW4: Mixed effects linear models


**Due**: Monday Nov. 22nd before class (printed copy at the class or earlier in my mailbox).

**Notes**:
  - For all the problems use at least 3500 iterations and 500 for burn-in.
  - You can get most of the code needed for the HW from Examples 1, 2 and 3 of this [entry](https://github.com/gdlc/STT465/blob/master/mixedEffects.md). Some additional modifications are needed, but the core of the code needed is there.

**Problem 1**. Compare goodness of fit and prediction accuracy of OLS and Bayesian regression.

Compute the correlation between predictions and observations in training and testing data sets using 30 training-testing partitions. Report: 
 1. A plot with training correlation for OLS versus training correlation for the Bayesian model and a 45-degree line added, 
 2. The same plot with testing correlation,
 3. A table with mean training and testing correlation per method and the proportion of times (over training-testing partitions) where the Bayesian gave higher correlation in each of the sets, 
 4. A summary of your findings (1 paragraph).

**Problem 2** Comparison of models based on cage versus models based on cage and molecular markers.

Repeat the analysis done in Problem 1 taking as the two models being compared a Bayesian model with cage as predictors and one with cage and SNPs as predictors.
