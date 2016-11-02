
## Example 1. OLS regression using 'cell means'

```R
 rm(list=ls())
 library(BGLR)
 data(mice)
 
 counts=table(mice.pheno$cage)
 tmp=names(counts)[counts>1]
 Y=mice.pheno[mice.pheno$cage%in%tmp,]
 y=Y$Obesity.BMI
 cage=factor(Y$cage)
 
 # Note 1: OLS on a set of orthogonal dummy variables (columns of Z) produces 'cell means'
  cageMeans=tapply(FUN=mean,X=y,INDEX=cage)
  fm=lm(y~cage-1)
  plot(cageMeans,coef(fm))
  
 # Note 2: SEs of OLS are inversely proportional to counts
  SE.lm=summary(fm)$coef[,2]
  n=counts[match(names(SE.lm),paste0('cage',names(counts)))]
  plot(n,SE.lm)
```

When the number of records per group (cage in our case) is small, the sampling variance of estimates is large and over-fitting becomes a risk. In these cases Bayesian estimators that shrink estimates towards a prior mean can aid. 

Before, we have assigned IID normal priors to effects with null mean and very large variance to avoid influence of the prior on inferences.
However, in this case we want the prior to influence estimates. Using a smaller prior variance can induce shrinkage of estimates towards zero. But what value should we choose for the prior variance? It turns out we can estimate the variance parameter from the data by simply treating this variance as unknown. For simplicity we will assign to it a scaled-inverse chi-square prior.


## Bayesian Shrinkage Estimation (Bayesian 'Ridge Regression')

