########################################################################
# A simple Gibbs Sampler for the mean and variance of a normal model   #
########################################################################

```R

## Toy simulation to test the program
 n=100
 y=rnorm(mean=125,sd=1,n=n)

## Hyper-parameters
 df0=0
 S0=0  # these two give the prior p(var)=1/var
 mu0=0
 varMu=1e6 # this gives flat prior for the mean
 
## Parameters that control the algorithm
 nIter=12000
 burnIn=2000

 
## Objects that will store samples
 mu=rep(NA,nIter)
 varE=rep(NA,nIter)
 
## Initializing
 varE[1]=var(y)
 mu[1]=mean(y)
 
## Objects that we will use
 Sy=sum(y)
 
## Sampler

 for(i in 2:nIter){
 	# sample the mean
 	 rhs=Sy/varE[i-1]+1/varMu

 	 C=n/varE[i-1]+1/varMu
 	 sol=rhs/C
 	 mu[i]=rnorm(sd=sqrt(1/C),mean=sol,n=1)
 	
 	# sampling the variance
 	 SS=sum((y-mu[i])^2)+S0
 	 DF=n+df0
 	 varE[i]=SS/rchisq(n=1,df=DF)
 	 print(i)
 }

 plot(varE)
 abline(h=var(y),lwd=2,col=2)
 plot(hist(varE,30))
 abline(v=var(y),lwd=2,col=2)
 plot(mu,varE)
```
