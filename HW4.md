### (1) Post-Gibbs Analyses

Using the gout data set fit a linear regression model with sex, age and race as effects and serum urate
as response (note the code for this is available in out GitHub repository, gibbsMLR.md). Collect a total of 20,000 samples.

   1.1. For each parameter in the model provide a trace plot and decide upon burn-in.
   Carry out the analyses requested below only based on samples that were not discarded as burn-in.
   
   
   1.2. Provide, for each parameter an estimate of the posterior mean, posterior SE and Monte-Carlo error
   (Hint: try summary on an mcmc object, created using the coda package).
   
   
   1.3. Provide, for each parameter, outo-correlation plots from lags 0 to lag 100.

   1.4. For each parameter provide an estimate of the number of independent samples that you have collected (hint: try effectiveSize())

   1.5. For each parameter in the model provide an estimate of the posterior density (try plot(density())). 
   Superimpose in each plot vertical lines indicating the estimated posterior mean and estimated 95% high-posterior density intervals (hint try HPDintervals())


### (2) Effects of shrinkage on prediction accuracy in High-dimensional regressions

For this data set we are going to use the wheat data set included in the BGLR package. 

After installing BGLR, you can access the data as follows


```r
 library(BGLR)
 data(wheat)
 X=scale(wheat.X) # genetic markers
 y=wheat.Y[,1]
 
 ## To make it faster I will take the first 500 columns
 
 X=X[,1:500]
 X=cbind(1,X) # adding a vector of 1's for the intercept

```

For this question the goal is to produce plots with prediction squared-correlation in the vertical axis versus
the value of the prior variance of effects in the horizontal axis for a training and a testing set. Use the following
code to partition your data into training and testing sets.


```r
 set.seed(195021)
 tst=sample(1:599,size=100,replace=F) 
 
 XTST=X[tst,]
 yTST=y[tst]
 
 XTRN=X[-tst,]
 yTRN=y[-tst]
```

For each of the values of varB in the following sequence `varB=exp(seq(from=log(.0001),to=log(10000),length=10))`:

    - Run a Gibbs sampler for 5500 iterations
    - Discard the frist 500 for burn-in
    - Using the remaining 5000, estimate the posterior mean of effects
    - Compute predictions in the TRN and in the TST set.
    

Report plots (one for the training and one for the testing set) with squared-correlation between predictions
and data in the vertical axis and values of varB in the horizontal axis.


Summarize your findings.
