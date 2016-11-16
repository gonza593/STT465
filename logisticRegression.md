# Logistic regression

Logistic regression is one of the most commonly used regression methods used for binary outcome. The regression function is defined
at the level of the logit

          
            log(theta_i/(1-theta_i))]xi'beta
            
 
 The following function evaluates the negative log-likelihood for a logistic regression.
 
 ```R
 
negLogLikLogistic=function(y,X,beta){
       eta=X%*%beta
       theta=exp(eta)/(1+exp(eta))
       
       out=ifelse(y==1,log(theta),log(1-theta))
       return(-sum(out))
}
 
 ```
 
 Now we show how the above function can be used in combination with `optim` to obtain Maximum Likelihood Estimates.
 
 We compare our MLE estimates with those obtained using the `glm` function.
 
 ```R
 load('~/GitHub/EPI853B/gout.RData')

 y=ifelse(Y$Gout=='Y',1,0)
 X=as.matrix(model.matrix(~UricAcid+Race+Sex+Age,data=Y))
 
 # centering all columsn of X does not change the ML estimates but facilitates convergence. 
  for(i in 2:ncol(X)){ X[,i]=X[,i]-mean(X[,i])}

 INT=log(mean(y)/(1-mean(y)))
 initialValues=c(INT, rep(0,ncol(X)-1))
 
 # Optim is a general purpouse optimization (minimization) function.
  fm=optim(fn=negLogLikLogistic,y=y,X=X,par=initialValues)
  names(fm$par)=colnames(X)
  
 # GLM performs ML estimation for logistic regression
  fm2=glm(y~X-1,family=binomial)
  
 # Comparison of results
  cbind(fm$par,coef(fm2))
 
 ```
