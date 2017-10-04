### Normal Model

We will discuss liklelihood and Bayesian inference for the normal model y~N(mu,V). We will begin with the case where the variance is known, 
then relax that assumption and discuss the model for the mean and variance. Then we will extend the model to the case where the mean
is parametrically modeled using a simple linear regression, so that yi~N(mui,V) with mui=mu+xib. Finally we will extend this to the case
where the mean is modeled using a multiple linear regression.

The first 2 cases are discussed in the book in chapter 5. Here I reproduce some of the scripts we develop in class. This is by no means
a substitute for the notes you take in class and the book.

#### (1) Mean with known variance

In the following example we examine the ML estimator and the mean and variance of the posterior density for data with the same mean
and three different sample sizes.

```r
## Data (same mean, varying sample size)
 yBar=15
 n=c(10,20,50)
 Vy=10 # variance of the data

## prior
 mu0=4 # prior mean
 Vmu=5  # prior variance

 ## MLE
  varMLE=Vy/n
  # Approx. 95% CI for the the Max. Likelihod estimator.
  CI=rbind(15+c(-1.96,1.96)*sqrt(varMLE[1]),
	   15+c(-1.96,1.96)*sqrt(varMLE[2]),
 	   15+c(-1.96,1.96)*sqrt(varMLE[3])
 	)
  show(CI)
  
 ## Posterior density
  tau2_data=1/Vy  # presicion of a single data-point
  tau2_mu=1/Vmu   # prior precision
 
  postMean=c(  (n[1]*yBar*tau2_data +tau2_mu*mu0)/( n[1]*tau2_data+tau2_mu) ,
 			  (n[2]*yBar*tau2_data +tau2_mu*mu0)/( n[2]*tau2_data+tau2_mu),
 			  (n[3]*yBar*tau2_data +tau2_mu*mu0)/( n[3]*tau2_data+tau2_mu)
 			)
 			
 			
   postVar=1/c(	n[1]*tau2_data+tau2_mu ,
 			n[2]*tau2_data+tau2_mu  ,
 			n[3]*tau2_data+tau2_mu
 			)
 			
    show(postMean)
    
    ## Posterior 95% credictility regions
    CR=cbind(qnorm(mean=postMean,sd=sqrt(postVar),p=.025),
          qnorm(mean=postMean,sd=sqrt(postVar),p=.975))
    show(CR)
 ```
 	
