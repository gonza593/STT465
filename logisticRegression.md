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
     theta=exp(
  }
 
 ```
