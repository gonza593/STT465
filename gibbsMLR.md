## Reading data to test the algorithm

```r
DATA=read.table('~/Desktop/gout.txt',header=F)
colnames(DATA)=c('sex','race','age','serum_urate','gout')

DATA$sex=factor(DATA$sex,levels=c('M','F'))
DATA$race=as.factor(DATA$race) 
```

#### Preparing incidence matrix 

```r
dF=ifelse(DATA$sex=='F',1,0) # a dummy variable for female
dW=ifelse(DATA$race=='W',1,0) # a dummy variable for male
 
# Incidence matrix for intercept and effects of sex, race and age
 X=cbind(1,dF,dW,DATA$age) 
 head(X)
 y=DATA$serum_urate
```

#### Gibbs

I wrapped up our Gibbs sampler into a function

```r
gibbsMLR=function(y,X,nIter=10000,varB,verbose=500){

 ## Hyper-parameters (these could be arguments to the function, for now we specify it here
  b0=0
  df0=4
  S0=var(y)*0.8*(df0-2)

 ## Objects to store samples
  p=ncol(X); n=nrow(X)
  B=matrix(nrow=nIter,ncol=p,0)
  varE=rep(NA,nIter)

 ## Initialize
  B[1,]=0
  B[1,1]=mean(y)
  b=B[1,]
  varE[1]=var(y)
  resid=y-B[1,1]
 
 ## Centering
  for(i in 2:ncol(X)){  X[,i]=X[,i]-mean(X[,i]) }

 ## Computing sum x'x for each column
  SSx=colSums(X^2)

  for(i in 2:nIter){
    # Sampling regression coefficients
    for(j in 1:p){
      C=SSx[j]/varE[i-1]+1/varB
      yStar= resid+X[,j]*b[j]  #y-X[,-j]%*%b[-j]
      rhs=sum(X[,j]*yStar)/varE[i-1]  + b0/varB
      condMean=rhs/C
      condVar=1/C
      b[j]=rnorm(n=1,mean=condMean,sd=sqrt(condVar))
      B[i,j]=b[j]  
      resid=yStar-X[,j]*b[j]
    }
    # Sampling the error variance  
    RSS=sum(resid^2)
    DF=n+df0
    S=RSS+S0
    varE[i]=S/rchisq(df=DF,n=1)

    if(i%%verbose==0){ cat(' Iter ', i, '\n') }
  }
 
  out=list(effects=B,varE=varE)
  return(out)
 }
```

**Use:**

```r
 samples=gibbsMLR(y,X,nIter=100,varB=10000)
 head(samples$varE) # samples for the error variance
 head(samples$effects)    # samples for the effects

```

### Effect of shrinkage on prediction accuracy in high-dimensional regressions

In this example we fit a linear model for wehat yield (y), using 1279 genetic markers (X) as covaraites. Sample size is 599. Therefore, in this problem the number of uknown regression coefficients exceeds sample size. The OLS estimate is not unique and the statistical properties of any of the soultions to OLS that can be obtained using a generalized inverse have very poor statistical properties. However, the Bayesian approach will enable us to get better estimates and even more a model that will have a reasonably good ability to predict grain yield from genetic markers in testing data sets.


In this example we split the data set into a training and a testing set, fit a Bayesian model to data in the training set and used the estimated effects to predict grain yield in the testing data set. To examine the role of shrinkage we fit the Bayesian model over a grid of values of the variance of effects, from very large (this gives a relatively un-informative prior and hence estimates that close to OLS solutions) to very small values of the variance (this will render estimates that strongly shrunked towards the prior mean (zero in the example). The example shows how the prior variance of effects can be used to induce shrinkage of estimates and how shrinkage can be used to arrive to models that will have good predictive performance.


The code can be run with our Gibbs sampler (the one printed above) or using BGLR. By default I run it with BGLR because the package is optimized and it is considerably faster than the Gibbs sampler that we have been using in class.


#### Fitting the models
```r

library(BGLR)
data(wheat)
X=wheat.X #  A matrix with genotypes of 599 lines of wheat (in rows) @ 1279 DNA markers (columns)
y=wheat.Y[,1] # grain yield of each of the lines


X=scale(wheat.X)/sqrt(ncol(X))

set.seed(195021)
tst=sample(1:599,size=150) # A testing set

XTST=X[tst,];XTST=cbind(1,XTST)
yTST=y[tst]

XTRN=X[-tst,];XTRN=cbind(1,XTRN)
yTRN=y[-tst]

varB=c(1000,100,10,5,2,1,.5,.1,.05,.01,.0001,1e-6) # Very large varB=>little shrinkage (in the limit Bayes=ML=OLS)

COR=rep(NA,length(varB))
BHat=matrix(nrow=ncol(X),ncol=length(varB))
varE=rep(NA,length(varB))

nIter=12000
burnIn=2000

for(i in 1:length(varB)){
	
	if(FALSE){ #done with our Gibbs sampler
		samples=gibbsMLR(y=yTRN,X=XTRN,varB=varB[i],nIter=nIter,verbose=100)
		bHat=colMeans(samples$effects[-(1:burnIn),])
		varE[i]=mean(samples$varE[-(1:burnIn)])
		BHat[,i]=bHat[-1]
	}else{ # done with BGLR that has a few additional optimizations and parts written in C
		DF=1e5 
		fm=BGLR(y=yTRN,ETA=list(list(X=XTRN[,-1],model='BRR',df0=DF,S0=varB[i]*(DF+2))),
          		nIter=nIter,burnIn=burnIn,verbose=F)
        BHat[,i]=fm$ETA[[1]]$b
        varE[i]=fm$varE
    	print(c(varB[i],fm$ETA[[1]]$varB))    
	
	}
	COR[i]=cor(yTST,XTST[,-1]%*%BHat[,i])
	print(COR[i])
}

# Examining estimates of effects

plot(numeric()~numeric(),xlab='varB',ylab='Estimated posterior mean',
     ylim=range(BHat[,1]),xlim=range(varB))

for(i in 1:1279){
	points(x=varB,y=BHat[i,],cex=.1,col=2)
	lines(x=varB,y=BHat[i,],lwd=.1,col=8)
}
```

#### Plots

```r
#  A plot in the -log scale for varB


plot(numeric()~numeric(),xlab='-log(varB)',ylab='Estimated posterior mean'
	,ylim=range(BHat[,1]),xlim=range(-log(varB)))

for(i in 1:1279){
	points(x=-log(varB),y=BHat[i,],cex=.1,col=2)
	lines(x=-log(varB),y=BHat[i,],lwd=.1,col=8)
}


plot(COR,x=log(varB),type='o',col=4)


```
