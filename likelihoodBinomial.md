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

Let's repeat the excercise many times, simulating each time N samples from conditonal IID `Bernoulli` and plot the likelihood function of each sample. To have a nicer plot, in each Monte-Carlo sample I scale the likelihood to have a max equal to one (this of course is just a scaling factor and therefore does not change the maximum or the shape of the likelihood).

```r
plot(numeric()~numeric(),xlim=c(0,1),ylim=c(0,1),ylab='Likelihood',xlab='theta')

# evaluating the log-lik for 500 samples (each of a size of 100) 
for(i in 1:500){
   Y=rbinom(size=1,n=100,prob=0.1)
   lines(x=x,y=exp(myLogLik(Y,x))/1e18,col='grey') # I scaled the likelihood
}
```
## A Monte Carlo Experiment to assess whether the Normal distribution that we assume (under CLT) is close to the actual distribution of the estimator over conceptual repeated sampling.


```r
nMC=1000
thetaHat=rep(NA,nMC)

theta=.5
nSamples=200

for(i in 1:nMC){
  Y=rbinom(prob=theta,size=1,n=nSamples)
  thetaHat[i]=mean(Y)
  print(i)
}


```
