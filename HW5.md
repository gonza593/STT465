 
## HW5: Bayesian Logistic Regression (via Metropolis-Hastings Algorithm)

Due Wednesday Dec. 7 in the class.

The data for this HW can be downloaded from the following [link](https://github.com/gdlc/EPI853B/blob/master/gout.RData).

The following example illustrates how to fit a logistic regression for Gout on Uric Acid via Maximum Likelihood.

 
```R
 load('~/GitHub/EPI853B/gout.RData') 

 y=ifelse(Y$Gout=='Y',1,0)
 x=Y$UricAcid
 
 # GLM performs ML estimation for logistic regression
 fm2=glm(y~x,family=binomial)
 summary(fm2)
```

Implement a Metropolis Hastings Algorithm for a Bayesian Logistic Regression. Use your algorithm to draw 15000 samples of the intercept and of the regression coefficient.

**Report**:
    - Trace plot for each of the parameters (how many iterations are you going to discard for burn-in?).
    - Estimated posterior means, SD and 95% credibility regions. 
    - Auto-correlation plots and estimated number of effective samples.
    - Are your estimates similar to ML?
    - Are you satisfied with your algorithm? If not, how could you improve it?
    
    
