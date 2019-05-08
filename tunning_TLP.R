###Tunning parameter selection of lambda and tau
######using five-fold cross validation
################
#source("vc_TLP.R")
#source("gTLP_2.R")
tunning_TLP<-function(x,tvar,y,N,p)
{
n<-length(y)
d<-dim(x)[[2]]
k<-ceiling(n/5)
h=n^(1/(2*p+3))

#r_tmp<-vc(x,tvar,y,x,tvar,N,p)
#alphahat<-r_tmp[[1]]
#tmp<-apply(alphahat^2,2,mean)
#tmp<-sqrt(tmp)

k1=20
k2=20
lgrid<-seq(0.01*log(n)*(n*h)^{-1/2},1.5*log(n)*(n*h)^{-1/2},length=k1) 
#tgrid<-c(0.01*log(n)*(n*h)^{-1/2},0.02*log(n)*(n*h)^{-1/2}) 
tgrid<-seq(0.01*log(n)*(n*h)^{-1/2},1.5*log(n)*(n*h)^{-1/2},length=k2) 

BIC<-matrix(0,k1,k2)

temp<-sample(n)
for(i in 1:k1)
{
for(j in 1:k2)
{
lambda<-lgrid[i]
tau<-tgrid[j]
result<-vc_TLP(x,tvar,y,x,tvar,N,p,lambda,tau)
betae<-result[[4]]
parp=sum(betae!=0)
yhat<-result[[3]]
BIC[i,j]=log(mean((y-yhat)^2))+log(n)*parp/n;
}
}
indx<-which(BIC==min(BIC),arr.ind=T)
lambda=lgrid[indx[1,1]]
tau=tgrid[indx[1,2]]
return(list(lambda,tau,lgrid,tgrid,BIC))
}