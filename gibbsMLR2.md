### We extend gibbsMLR() by addin an step to sample varB.

```r

gibbsMLR2=function(y,X,nIter=10000,verbose=500){

 ## Hyper-parameters
  b0=0
  df0=4
  S0=var(y)*0.8*(df0-2)
  
  df0b=4
  S0b=0.1

 ## Objects to store samples
  p=ncol(X); n=nrow(X)
  B=matrix(nrow=nIter,ncol=p,0)
  varE=rep(NA,nIter)
  varB=rep(NA,nIter)
  
 ## Initialize
  B[1,]=0
  B[1,1]=mean(y)
  b=B[1,]
  varE[1]=var(y)
  varB[1]=S0b/(df0b+2)
  resid=y-B[1,1]
 print(varB[1])

 ## Centering
  for(i in 2:ncol(X)){  X[,i]=X[,i]-mean(X[,i]) }

 ## Computing sum x'x for each column
  SSx=colSums(X^2)

  for(i in 2:nIter){
    # Sampling regression coefficients
    for(j in 1:p){
      C=SSx[j]/varE[i-1]+1/varB[i-1]
      tmp=sum(X[,j]*resid)+SSx[j]*b[j]  #y-X[,-j]%*%b[-j]
      rhs=tmp/varE[i-1]  + b0/varB[i-1]
      condMean=rhs/C
      condVar=1/C
      b_old=b[j]
      b[j]=rnorm(n=1,mean=condMean,sd=sqrt(condVar))
      B[i,j]=b[j]
      resid=resid-X[,j]*(b[j]-b_old) # updating residuals
    }

    # Sampling the error variance
    RSS=sum(resid^2)
    DF=n+df0
    S=RSS+S0
    varE[i]=S/rchisq(df=DF,n=1)

    # Sampling the variance of effects
    DF=p+df0b
    S=S0b+sum((b-b0)^2)
    varB[i]=S/rchisq(n=1,df=DF)

    if(i%%verbose==0){ cat(' Iter ', i, '\n') }
  }
  out=list(effects=B,varE=varE,varB=varB)
  return(out)
 }
