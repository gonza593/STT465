DATA=read.table('~/GitHub/STT465/gout.txt',header=F)
colnames(DATA)=c('sex','race','age','serum_urate','gout')

DATA$sex=factor(DATA$sex,levels=c('M','F'))
DATA$race=as.factor(DATA$race) 

dF=ifelse(DATA$sex=='F',1,0) # a dummy variable for female
dW=ifelse(DATA$race=='W',1,0) # a dummy variable for male
 
# Incidence matrix for intercept and effects of sex, race and age
 X=cbind(1,dF,dW,DATA$age) 
 head(X)
 y=DATA$serum_urate

## Hyper-parameters
 varB=1000
 b0=0
 varE=1.8

## Number of iterations
 nIter=10000

## Objects to store samples
 p=ncol(X)
 B=matrix(nrow=nIter,ncol=p,0)

 SSx=colSums(X^2)

## Initialize
 B[1,]=0
 B[1,1]=mean(y)
 b=B[1,]

 for(i in 2:nIter){

  for(j in 1:p){

    C=SSx[j]/varE+1/varB

    yStar=y-X[,-j]%*%b[-j]

    rhs=sum(X[,j]*yStar)/varE  + b0/varB
    
    condMean=rhs/C
    condVar=1/C

    b[j]=rnorm(n=1,mean=condMean,sd=sqrt(condVar))
    B[i,j]=b[j]  
}
 print(i)
 }


