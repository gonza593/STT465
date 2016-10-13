### HW 3:  Linear Regression
**Due Friday Oct. 21st** (leave a printed copy of your HW in my mailbox in the STT department, Wells C, 4th floor).


1. For a simple linear regression model y=mu+xb+e, write down:  
   - (1.a) The sampling model assuming IID normal errors
   - (1.b) The prior distribution of the intercept, regression coefficient and variance assuming normal independent priors for intercept and regression coefficient and scaled-inverse chi-square for the variance.
   
2. Suppose that a-priori you believe that the error variance should be around ~2, if the DF of parameter is equal to 4, what should the scale parameter be so that the 
  - (2.a) the prior mean equals 2
  - (2.b) the prior mode equals 2
3. Suggest values for the prior mean and variance of the priors of the intercept and of the regression coefficient that will yield negligable influence of the prior on inferences.

4. Derive the fully conditional distributions of each of the parameters of the model (show your derivations step by step)

5. Implement a Gibbs Sampler using the fully conditionals you derived

6. Using the [gout](https://github.com/gdlc/STT465/blob/master/gout.txt) data set run 1200 iterations of your sampler for a regression of LDL cholesterol levels on sex. Discard the first 200 iterations for burn in.
  - (6.a) Report trace plots for each of the parameters
  - (6.b) Report density plots for each of the parameters
  - (6.c) What is the estimated difference in LDL cholesterol between male and female? Report an estimate and a 95% posterior credibility region.
    
    
Hint: you can introduce sex as a regressor by creating a dummy variable that takes 1 for one of the genders (e.g., female) and zero otherwise.
 
 

