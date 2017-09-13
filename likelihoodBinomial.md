### The likelihood portion of the Beta-Binomial model.


Consider a sample of N (conditional) IID bernouly random variables, Yi={0/1}, i=1,...,N. 

**Evaluating the likelihood function**

```r
  theta=0.2 # uknown true succses probability
  nSamples=20 # sample size
  
  Y=rbinom(size=1,n=nSamples,prob=theta)
  head(Y)
  
  mean(Y) # the ML estimator
  
  myLogLik=function(Y,theta){
     N=length(Y)
     N1=sum(Y)
     N2=N-N1
     logLik=N1*log(theta)+N2*log(1-theta)
     return(logLik)
  }
  
  x=seq(from=1/100,to=99/100,length=100 ) # a grid of values of theta
  y=exp(myLogLik(Y,x))
  plot(y~x,type='l',ylab='Likelihood', xlab='Success Probability'); abline(v=mean(Y))
  
  ## Another option
  y2=dbinom(x=sum(Y),size=nSamples,prob=x) # note there is a proportionality constant difference between y1 and y2, but the ML will be the seame
```

## The likelihood as a random function

Let's repeat the excercise many times, simulating each time N samples from conditonal IID `Bernoulli` and plot the likelihood function of each sample.

```r


```
