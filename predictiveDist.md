## Computing the predictive distribution using MC methods

**Example 1: Binomial**

Here becasue the future data (yF) is a dummy variable, the predictive distribution must be Bernoulli, we just
need to find the success probability. It turns out that the success probability (see derivation in class or in the book, p. 40) is the 
posterior mean of theta. The example below compares the predictive distribution derived using this result with one obtained using MC methods.


```R
 B=5e6
 
 shape1=1.5
 shape2=1.5
 N=20
 xBar=.2

 # posterior density is Beta with these shape paramerters
  post.a=N*xBar+shape1
  post.b=N*(1-xBar)+shape2
 
 theta=rbeta(shape1=post.a,shape2=post.b,n=B)
 yF=rbinom(p=theta,n=B,size=1)
 
 mean(yF) # estimated success probability obtained by averaging over possible realizations of theta
 (N*xBar+shape1)/(N+shape1+shape2) # analythical solution

```

**Example 2: Poisson-Gamma:**

In this case the predictive density is Negative Binomial, but let's try to obtain it using MC methods.


```R
 # Simulating data from Poisson
  N=20
  x=rpois(lambda=4,n=N)
 
 # Setting the prior, say that our prior expectation is E[lambda]=3 and CV[lambda]=1 (large variance relative to mean)
  prior.mean=3
  prior.CV=1
  
  prior.shape=(1/prior.CV)^2
  prior.rate=prior.shape/prior.mean

 # The posterior distribution is gamma with these parameters
  post.shape=prior.shape+sum(x)
  post.rate=N+prior.rate
 
 # Let's generate the predictive distribution
 
  B=1e6
  samplesLambda=rgamma(n=B,shape=post.shape,rate=post.rate)
  yF=rpois(n=B,lambda=samplesLambda)
  
  var(rpois(n=B,lambda=post.shape/post.rate))

```
