# Data
 trueMu=10
 trueVar=1.5
 n=23
 y=rnorm(n=n,mean=trueMu,sd=sqrt(trueVar))

# Prior hyper-parameters
 v=3
 S=6  # will discuss later on how to chose this 
 mu0=0
 varMu=10000 # flat prior

# Parameters of the sampler
 nIter=10000

# Storage for the samples
 mu=rep(NA,nIter)
 varY=rep(NA,nIter)
 
# Initialize
 meanY=mean(y)

 mu[1]=meanY
 
 DF=v+n
 RSS=sum( (y-mu[1])^2)
 tmpS= S + RSS
 varY[1]=tmpS/rchisq(n=1,df=DF)# Initialize

# Gibbs Sampler
 for(i in 2:nIter  ){
    
    # Sampling mu given data and variance
     condVar=1/( n/varY[i-1] + 1/varMu )
     w=(n/varY[i-1])*condVar
     condMean=w*meanY+(1-w)*mu0
     mu[i]=rnorm(n=1,mean=condMean,sd=sqrt(condVar))
    
    # sampling variance given data and mu
     RSS=sum( (y-mu[i])^2)
     tmpS= S + RSS
     varY[i]=tmpS/rchisq(n=1,df=DF)# Initialize
    print(i)
 }



