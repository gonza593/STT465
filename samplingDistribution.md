## Simple illustration of the (Frequentist) concept of the sampling distribution of an estimator.

Suppose we want to estimate the mean height in a population. Let's assume that height follows a Normal distribution with mean *mu* and variance
*V*. This distribution is of cours uknown to us. To estimate the **population mean** this is our plan: we will collect N samples from the population, 
compute the **sample mean** and use this as our **estimate**. The code below illusrates this.


```r
 n=10
 mu=175
 V=40
 mySample=rnorm(mean=mu,sd=sqrt(V),n=n)
 sum(mySample)/length(mySample)
 mean(mySample)

```

What would happen if we go again to the population and collect another sample?



```r
 n=10
 mu=175
 V=40
 mySample=rnorm(mean=mu,sd=sqrt(V),n=n)
 sum(mySample)/length(mySample)
 mean(mySample)

```

The estimated value varies over conceptual repeated sampling....


**The Sampling Distribution**

Consider now sampling an infinte number of times....


```r
   n=10
   mu=175
   V=40
 
  nRep=1000 # this is not infinite but it is large enough...
  estimates=rep(NA,nRep)
  for(i in 1:nRep){
    mySample=rnorm(mean=mu,sd=sqrt(V),n=n)
    estimates[i]=mean(mySample)
  }
  hist(estimates,30)
```

[return](https://github.com/gdlc/STT465/)
