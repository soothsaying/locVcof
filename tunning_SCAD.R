###Tunning parameter selection of lambda and tau
######using five-fold cross validation
################
source("vc_SCAD.R")
tunning_SCAD<-function(x,tvar,y,N,p)
{
n<-length(y)
d<-dim(x)[[2]]
k<-ceiling(n/5)
h=n^(1/(2*p+3))
k1=20
lgrid<-seq(0.01*log(n)*(n*h)^{-1/2},1.5*log(n)*(n*h)^{-1/2},length=k1) 
BIC<-rep(0,k1)
a=3.7

temp<-sample(n)
for(i in 1:k1)
{
lambda=lgrid[i]
result<-vc_SCAD(x,tvar,y,x,tvar,N,p,lambda,a)
betae<-result[[4]]
parp=sum(betae!=0)
yhat<-result[[3]]
BIC[i]=log(mean((y-yhat)^2))+log(n)*parp/n;
}

ind=which.min(BIC)
return(lgrid[ind])
}
