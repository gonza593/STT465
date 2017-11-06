


DATA$sex=factor(DATA$sex,levels=c('M','F'))
fm=lm(serum_urate~sex,data=DATA)
summary(fm)
fm$coef
fm$coef[2]
fm$coef[2]/0.02769
X=cbind(1,ifelse(DATA$sex=='F'))
X=cbind(1,ifelse(DATA$sex=='F',1,0))
head(X)
y=DATA$serum_urate
Xy=t(X)%*%y
XtX=t(X)%*%X
solve(XtX,Xy
)
fm$coef
bHat=solve(XtX,Xy)
eHat=y-X%*%bHat
var(eHat)
sum(eHat^2)/(nrow(DATA)-2)
summary(fm)
1.36^2
varEHat=sum(eHat^2)/(nrow(DATA)-2)
XXInv=solve(XtX)
diag(XXInv)*varEHat
sqrt(diag(XXInv)*varEHat)
