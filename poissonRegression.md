# Poisson Regression: Likelihood and Bayesian Analyses

### Loading the data into R

```r
 DATA=read.table('~/Desktop/crab.txt',header=T)
 DATA$color=factor(DATA$color)
 DATA$spine=factor(DATA$spine)
```

### Maximum Likelihood estimation using glm
```r
 fm0=glm(nSatellites~color+spine+width,family=poisson(link=log),data=DATA)
```
 
## A function to evaluate the negative log-likelihood
```r
negLogLikPoisson<-function(y,X,b){
 	eta<-X%*%b
 	lambda=exp(eta)
 	logLik=dpois(x=y,lambda=lambda,log=T)
 	return(-sum(logLik))
 }
```

## Fitting the model using `optim`
```r
 y=DATA$nSatellites
 X=as.matrix(model.matrix(~color+spine+width,data=DATA))[,-1]
 X=scale(X,center=T,scale=F)
 X=cbind(1,X)
 init=c(log(mean(y)),rep(0,ncol(X)-1)) #
 fm=optim(fn=negLogLikPoisson, y=y,X=X,par=init)
 fm$convergence # it did not converge

 #using previous fit as initial values and increasing the maximum number of iterations
 fm=optim(fn=negLogLikPoisson, y=y,X=X,par=fm$par,control=list(maxit=100000)) 
 fm$convergence
 
 round(cbind(fm0$coef,fm$par),5)
 
```


### Bayesian Analyses using Metropolis
 

**Setting the stage**
```r

 n=nrow(X)
 p=ncol(X)
 b0=rep(0,p)
 varB=rep(1000,p)
 V=.01 # variance of the distribution used to generate candidates
 
 logP<-function(y,X,b,varB,b0){
 	logLik<- -negLogLikPoisson(y=y,X=X,b=b)
 	logPrior<- sum(dnorm(x=b,mean=b0,sd=sqrt(varB),log=T))
 	return(logLik+logPrior)
 }
 
 nIter=60000
 burnIn=10000
 
 B=matrix(nrow=nIter,ncol=p,NA)
 B[1,]=0
 B[1,1]=log(mean(y))
 
 accept<-rep(NA,nIter)
 accept[1]=TRUE
 current<-B[1,]
 
```

**Sampler**

```r
 for(i in 1:nIter){
 	
 	candidate=rnorm(n=p,mean=current,sd=sqrt(V))
 	
 	r=min(1,exp(logP(y,X,candidate,varB,b0)-logP(y,X,current,varB,b0)))
 	accept[i]<- runif(1)<=r
 	
 	if(accept[i]){
 		current<-candidate
 	}
 	B[i,]=current
 	print(i)
 }
 # Comparing posterior means with Maximum Likelihood estimates
 cbind(fm$par,colMeans(B[-(1:10000),]))

```


**Tasks**:
  - Carry out diagnostics to determine: burnin period, Monte-Carlo error, auto-correlation of samples, effective sample size.
  - Once you decide on the burn-in period discard the samples within the burn-in period, and report estimated mean, estimated SD, estimated posterior credibility regions and estimated posterior density for each coefficient.
  - Interpret your results.
  - Summarize your findings.


