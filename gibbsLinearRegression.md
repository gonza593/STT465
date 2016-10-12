### A Gibbs Sampler for a linear regression model

```R

## Toy simulation to test the program
 n=100
 x=rnorm(n=n)
 int=100
 beta=2
 y=int+x*beta+rnorm(n)
 
 
  X=cbind(1,x)
  n=nrow(X)
  p=ncol(X)
 
## Hyper-parameters
 df0=0
 S0=0  # these two give the prior p(var)=1/var
 b0=rep(0,p)
 varB=rep(1e6 ,p)# this gives flat prior for the intercept and regression coef.
 
 
## Parameters that control the algorithm
 nIter=100
 
 
## Objects that will store samples
 B=matrix(nrow=nIter,ncol=p,NA)
 B[1,1]=mean(y) #initializing the intercept
 B[1,-1]=0 # now regression coef
 varE=rep(NA,nIter)
 varE[1]=var(y)
 SSx=colSums(X^2) # we will need this in the sampler
  
## Sampler

 for(i in 2:nIter){
 	b=B[i-1,]
  for(j in 1:p){ 
    yStar=y-X[,-j,drop=F]%*%b[-j]
    rhs=sum(X[,j]*yStar)/varE[i-1]+b0[j]/varB[j]
    C=SSx[j]/varE[i-1]+1/varB[j]
    sol=rhs/C
    b[j]=rnorm(n=1,mean=sol,sd=sqrt(1/C))
  }
  B[i,]=b
  # sampling the error variance
 	 error=y-X%*%b
 	# sampling the variance
 	 SS=sum(error^2)+S0
 	 DF=n+df0
 	 varE[i]=SS/rchisq(n=1,df=DF)
 	 print(i)
 }

 plot(varE)
 plot(B[,1])
 plot(B[,2])
 
```
