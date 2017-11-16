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

####Gibbs

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
  SSx=colSums(X^2)

 ## Initialize
  B[1,]=0
  B[1,1]=mean(y)
  b=B[1,]
  varE[1]=var(y)
  resid=y-B[1,1]
 
 ## Centering
  for(i in 2:ncol(X)){  X[,i]=X[,i]-mean(X[,i]) }
 
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
