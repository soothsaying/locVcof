###Tunning parameter selection of lambda and tau
######using five-fold cross validation
################
source("gLasso.R")
source("vc_LASSO.R")

tunning_AdLASSO<-function(x,tvar,y,N,p,weight)
{
n<-length(y)
d<-dim(x)[[2]]
k<-ceiling(n/5)
h=n^(1/(2*p+3))
#weight<-rep(1,d)

#r_tmp<-vc(x,tvar,y,x,tvar,N,p)
#alphahat<-r_tmp[[1]]
#tmp<-apply(alphahat^2,2,mean)
#tmp<-sqrt(tmp)

k1=20
k2=1
lgrid<-seq(0.01*(n*h)^{-1/2},0.1*(n*h)^{-1/2},length=k1) 
#tgrid<-c(0.01*log(n)*(n*h)^{-1/2},0.02*log(n)*(n*h)^{-1/2}) 
#tgrid<-seq(0.01*(n*h)^{-1/2},1.5*(n*h)^{-1/2},length=k2) 
CVE_LASSO<-rep(0,k1)

temp<-sample(n)
for(i in 1:k1)
{
lambda<-lgrid[i]
for(r in 1:5)
{
indtemp=((r-1)*k+1):(r*k)
indm<-temp[-indtemp]
indp<-temp[indtemp]
xm<-x[indm,]
tm<-tvar[indm]
ym<-y[indm]
xp<-x[indp,]
tp<-tvar[indp]
yp<-y[indp]
result<-vc_LASSO(xm,tm,ym,xp,tp,N,p,lambda,weight)
ypred_LASSO=result[[3]]
CVE_LASSO[i]<-CVE_LASSO[i]+mean((ypred_LASSO-yp)^2)
#print(result[[2]])
}
#print(lambda)
}

#print(CVE_TLP)

ind=which.min(CVE_LASSO)
return(lgrid[ind])
}
