### HW 5: Bayesian Analysis of the Logistic Regression

**Topics**:
  - Logistic Regression
  - Metropolis Algorithm
  - Analyses of Monte Carlo Samples
  
 **Due**: December 6th (Wed.) in class.
  
 **(1) Maximum Likelihood Estimaton and Interpretation of Parameter Estimates**
 
Using the [gout]() data set Fit a logistic regression with gout as response, sex, age, race and serum urate as predictors using the `glm()` function.


1.1. Report parameter estimates, SEs and p-values


1.2. Summarize your findings


1.3. Using the parameter estimates report the estimated probability of developing gout fo reach of the following cases



**(2) Bayesian Analysis**

Use the Gibbs sampler developed in class to fit the logistic regression. Collect 22,000 samples, discard the frist 2,000 for burn in.


2.1. Report parameter estimates, posterior standard deviation and 95% posterior credibility regions



2.2. Report for each coefficient the trace plot, and a table with numbrer of effective samples, and MC estandard error.


2.3. Summarize your finding


2.4. Use the samples collected to estimate the posterior distribution of the probability of developing gout for e each of these cases
(report one plot with the estimated posterior desnity of each of these probabilities).


|----------|----------|----------|----------|
| Sex    |  Race | Age |  Serum Urate |
|----------|----------|----------|----------|
| Sex    |  Race | Age |  Serum Urate |
|----------|----------|----------|----------|




