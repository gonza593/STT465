

### Metropolis Algorithm: a very simple example


Here we show how the Metropolis algorithm can be used to dra samples from a target distribution (`p()`) generating candidates 
from another distribution (q`(x)`). For the metropolis algorithm `q(x)` must be symmetric and have the same support than that of `p(x)`.

In the example below `x~N(0,1)`; we draw `n` samples from this distribution and use these samples for comparison with samples drawn
with a Metropolis algorithm.

```r
 n=100000    # number of samples we want to draw
 x=rnorm(n)  # samples from the target distribution (used for comparison only).
 
 z= rep(NA,n) 
 z[1]=-1
 for(i in 2:n){
 	candidate=runif(min=-5,max=5,n=1) # symmetric candidate generator
 	r=dnorm(candidate)/dnorm(z[i-1])
 	accept<-runif(1)<r
 	if(accept){
 		z[i]=candidate
 	}else{
 		z[i]=z[i-1]
 	}
 	print(i)
 	
 }
 
 plot(density(x),col=4)
 tmp=density(z)
 lines(x=tmp$x,y=tmp$y,col=2)
 
```


### Logistic Regression

Let `yi` be a Bernoulli variable (`yi` can take values 0/1) with success probability `p(Yi=1)=pi`.
In a logistic regression the succes probability ot each subject is modeled using a linear regression. Specifically,
`log(pi/(1-pi))=x1i*b1+x2i*b2+...`.


**Fitting a Logistic Regression using glm**

```r	
# Data
  DATA=read.table('~/Desktop/gout.txt',header=F)
  colnames(DATA)=c('sex','race','age','serum_urate','gout')
  DATA$gout=ifelse(DATA$gout=="Y",1,0)

# Incidence matrix for effects
  X=as.matrix(model.matrix(~sex+race+age+serum_urate,data=DATA))[,-1]
  X=scale(X,center=T,scale=F)
  
# glmm fit
  fm=glm(gout~sex+age+serum_urate,family=binomial(link=probit),data=DATA)
  summary(fm)

```
**A function that evaluates the log-likelihood,function.

```r
  negLogLik=function(y,X,b){
  	ETA=X%*%b
  	theta=exp(ETA)/(1+exp(ETA))
  	
  	logLik<-ifelse(y==1,log(theta),log(1-theta))
  	return(-sum(logLik))
  }
```  
  
**Fitting a logistic regression using `optim()`.**

```r

 
  X=cbind(1,X)
  yBar=mean(DATA$gout)
  init=c(log(yBar/(1-yBar)),rep(0,ncol(X)-1))  
  
  fm2=optim(X=X,y=DATA$gout,fn=negLogLik,par=init)

  
```



