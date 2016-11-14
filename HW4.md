# HW4: Mixed effecst linear models


**Due**: Monday Nov. 22nd before class (printed copy at my mailbox).


**Problem 1**. Using the mice data set fit a linear regression of BMI on an intercept plus cage effects. Througout the problem 
use 2200 iterations and a burn-in of 200. 


1.1. Fit the model to the entire data set and report for the error variance and the variance of cage effects: (i) trace and density plot 
(i) estimated posterior means, SD and 95% credibility regions, (iii) estimated effective number of samples.

1.2. Conduct a 5-fold cross-validation with assigment of individuals within cage to folds. Report for each fold the training and testing
correlation.


**Problem 2** Repeat the analysis of question 1 using both cages and SNPs as predictors. Estimate a variance specific to cage and one specific to SNPs.

2.1. Report the same items requested in 1.1 including in this case the variance of SNP effects as well.

2.2. Conduct the same analysis and report the same items requested in 1.2 for the model using both cage and SNPs as effect.s


**3** Summarize your conclusions. Does adding SNP improve the model? Why yes, why not?
