logLikPoiss=function(N,xBar,lambda){
    x1=N*xBar*log(lambda)
    x2=N*lambda
    return(x1-x2)
}


lambda=4
N=100

xBar=rnorm(sd=sqrt(lambda/N),mean=lambda,100)

lambda=sort(runif(1000)*20)

LOG_LIK=matrix(nrow=length(lambda),ncol=length(xBar))

for(i in 1:ncol(LOG_LIK)){
	LOG_LIK[,i]=logLikPoiss(N=10,xBar=xBar[i],lambda)
}

tmp=range(LOG_LIK)
plot(y=LOG_LIK[,1],x=lambda,type='l',ylim=tmp)
for(i in 1:ncol(LOG_LIK)){
   lines(x=lambda,y=LOG_LIK[,i])
}

