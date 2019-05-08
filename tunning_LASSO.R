###Tunning parameter selection of lambda and tau
######using five-fold cross validation
################
source("gLasso.R")
source("vc_LASSO.R")

tunning_LASSO<-function(x,tvar,y,N,p,weight)
{
n<-length(y)
d<-dim(x)[[2]]
k<-ceiling(n/5)
h=n^(1/(2*p+3))
weight<-rep(1,d)

#r_tmp<-vc(x,tvar,y,x,tvar,N,p)
#alphahat<-r_tmp[[1]]
#tmp<-apply(alphahat^2,2,mean)
#tmp<-sqrt(tmp)

k1=20
k2=1
lgrid<-seq(0.5*(n*h)^{-1/2},2*(n*h)^{-1/2},length=k1) 
#tgrid<-c(0.01*log(n)*(n*h)^{-1/2},0.02*log(n)*(n*h)^{-1/2}) 
#tgrid<-seq(0.01*(n*h)^{-1/2},1.5*(n*h)^{-1/2},length=k2) 
BIC<-rep(0,k1)

for(i in 1:k1)
{
lambda<-lgrid[i]
result<-vc_LASSO(x,tvar,y,x,tvar,N,p,lambda,weight)
betae<-result[[4]]
parp=sum(betae!=0)
yhat<-result[[3]]
BIC[i]=log(mean((y-yhat)^2))+log(n)*parp/n;
}


ind=which.min(BIC)
return(lgrid[ind])
}
