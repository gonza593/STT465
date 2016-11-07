
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

The variance of effects can be treated as uknown, in the same way we treated the error variance as uknown. The following sampler allows us to do that. The algorithm involves minimal changes relative to the one we used for the case of a model with a 'flat' prior for effects. I took the code of the [Gibbs Sampler for Multiple Linear Regression](https://github.com/gdlc/STT465/blob/master/gibbsLinearRegression.md) and made modifications so that we can group predictors into sets. Each set has it's own variance of effects. For some we can specify `type='fixed'` in which case the prior variance is assinged a large value and not updated, or `type='random'` in which case the variance parameter is sampled from the posterior density.

#### Gibbs Sampler

```R
GIBBS.MM=function(y,X,group,type,nIter,burnIn){
 n=nrow(X)
 p=ncol(X)
## Centering all columns except the 1st (intercept)
  for(i in 2:p){ X[,i]=X[,i]-mean(X[,i]) }
  SSx=colSums(X^2) # we will need the sum of squares in the sampler
  
    
## Hyper-parameters
 nGroups=length(unique(group))
 groupSize=table(group)
 df0.b=rep(0,nGroups)
 S0.b=rep(0,nGroups)  # these two give the prior p(var)=1/var
 b0=rep(0,nGroups)
 df0.e=0
 S0.e=0
 
## Objects that will store samples
 B=matrix(nrow=nIter,ncol=p,NA)
 B[1,1]=mean(y) #initializing the intercept
 B[1,-1]=0 # now regression coef
 error=y-mean(y) #*#
 varE=rep(NA,nIter)
 varE[1]=var(y)
 
 varB=matrix(nrow=nIter,ncol=nGroups)
 colVar=apply(FUN=var,X=X,MARGIN=2)
 
 Vy=var(y)
 Vmodel=.5*Vy
 Vb=Vmodel/sum(colVar)
 for(i in 1:nGroups){
 	if(type[i]=='fixed'){
 		varB[,i]=1e6
 	}else{
 	     varB[1,i]=Vb
 	}
 }

## Sampler

 for(i in 2:nIter){
  # Sampling effects
  b=B[i-1,]
  for(j in 1:p){ 
    error=error+X[,j]*b[j] #*#
     rhs=sum(X[,j]*error)/varE[i-1]+b0[group[j]]/varB[i-1,group[j]] #*#
     C=SSx[j]/varE[i-1]+1/varB[i-1,group[j]] #*#
     sol=rhs/C
     b[j]=rnorm(n=1,mean=sol,sd=sqrt(1/C))
    error=error-X[,j]*b[j]
  }
  B[i,]=b
  # sampling the error variance
 	 SS=sum(error^2)+S0.e
 	 DF=n+df0.e
 	 varE[i]=SS/rchisq(n=1,df=DF)
 	 print(i)
  # sampling variance of effects #*#
  for(j in 1:nGroups){
  	if(type[j]!='fixed'){ 
  		tmp=b[group==j]-b0[j]
  		SS=sum(tmp^2)+S0.b[j]
  		DF=groupSize[j]+df0.b
  		varB[i,j]=SS/rchisq(n=1,df=DF)
  	}
  }
 }
 
 out=list(varE=varE,B=B,varB=varB)
 return(out)
}
```

## Example



