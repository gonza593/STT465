## Bayesian Analysis of the Logistic Regression


Logistic regression is used to estimate the effects of a set of predictors (e.g., sex, age) on the success probability of a binary outcome (Yi). This is done
by linking the success probability, P(Yi=1|xi,b), with a linear predictor using a non-linear transformation  (the *logit*). 

  * *Success probability* `P(Yi=1)=pi`
  * *Linear predictor* `LPi=xi'b` where `xi` is a vector of covariantes and `b` is  a vector of effects
  * *Odds* `pi/(1-pi)`.
  * *Logit* the logaritm of the odds, `log(pi/(1-pi))=LPi=xi'b`.
  * Solving the Logit for `pi` we get `pi=exp(LPi)/(1+LPi)`.
 
 Using the terms defined above we can evaluate the likelihood functio of a logistic regression using the folowing
 
 
 **Log-likelihood of a Logistic Regression**
 
 ```r
  loglik=function(y,X,b){
     LP=X%*%b
     theta=exp(LP)/(1+exp(LP))
     logLik=ifelse(y==1,log(theta),log(1-theta))
     return(sum(logLik))
  }
 
 ```

**Fitting a logistic regression via maximum likelihood using `glm`**

```r
rm(list=ls())
 
DATA=read.table('~/Desktop/gout.txt',header=F)
colnames(DATA)=c('sex','race','age','serumUrate','gout')
DATA$gout=ifelse(DATA$gout=="Y",1,0)
 
### Logistic Regression
 
 
# Incidence matrix for effects
Z=as.matrix(model.matrix(~sex+age+race+serumUrate,data=DATA))[,-1]
Z=scale(Z,center=T,scale=F) # centering to improve mixing and convergence
 
# Maximum Likelihood using glm
fm=glm(gout~Z,data=DATA,family='binomial')
summary(fm)

```

**Solving the maximum likelihood problem using a general-purpouse optimization function (`optim()`)**

```r
# Maximum Likelihood using optim
  negLogLik=function(y,X,b){
                 Xb=X%*%b
                 theta=exp(Xb)/(1+exp(Xb))
                 logLik=sum( y*log(theta)+(1-y)*log(1-theta))
                 return(-logLik)
  }
 
  iniInt=log(mean(DATA$gout)/(1-mean(DATA$gout)))
 
  fm2=optim(fn=negLogLik,y=DATA$gout,X=cbind(1,Z),par=c(iniInt,rep(0,ncol(Z)))
                                                ,control=list(maxit=10000,reltol=1e-8)) # these parameters control the 
                                                                                        # algorithm try help(optim) 
  fm2$convergence # 0 means it converged
  cbind(fm$coef,fm2$par)

```

