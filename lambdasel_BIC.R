############cross validation function
lambdasel_BIC<-function(x,y,tvar,N,p,weight,flag)
{
#lambdav<-seq(0.0001,0.1,length=200) #LASSO

n<-length(y)
m<-ceiling(n/5)
lambdav<-seq(0.01*(n/5)^{-1/2},1.5*(n/5)^{-1/2},length=200) 

###TLP
if(flag==4){
tauv<-seq(0.001,2,length=50)
cvr<-matrix(0,nrow=length(lambdav),ncol=length(tauv))
for(i in 1:length(lambdav))
{
for(j in 1:length(tauv))
{
cvij=0
for(k in 1:5)
{
ind_pred<-indt[((k-1)*m+1):(k*m)]
xpred=x[ind_pred,]
xest=x[-ind_pred,]
yest=y[-ind_pred]
ypred=y[ind_pred]
ttest=tvar[-ind_pred]
tpred=tvar[ind_pred]
result<-splinevcf(as.matrix(xest),yest,ttest,N,p,tpred,xpred,lambdav[i],tauv[j],0,flag)
yest<-result[[3]]
cvij=cvij+mean((ypred-yest)^2)
print(mean((ypred-yest)^2))
}
cvr[i,j]=cvij
}
}
indx<-which(cvr==min(cvr),arr.ind=T)
lambda=lambdav[indx[1,1]]
tau=tauv[indx[1,2]]
return(list(lambda,tau,lambdav,tauv,cvr))
}

else
{
BIC<-rep(0,length(lambdav))
for(i in 1:length(lambdav))
{
result<-splinevcf(x,y,tvar,N,p,tvar,x,lambdav[i],0,weight,flag)
betae<-result[[1]]
parp=sum(betae!=0)
yhat<-result[[3]]
BIC[i]=log(mean((y-yhat)^2))+log(n)*parp/n;
}
indx<-which.min(BIC)
lambda=lambdav[indx]
return(list(lambda,lambdav,BIC))
}

}




############cross validation function
lambdasel_CV<-function(x,y,tvar,N,p,weight,flag)
{
n<-length(y)
m<-ceiling(n/5)
indt<-sample(n)
lambdav<-seq(0.01*(n/5)^{-1/2},1.5*(n/5)^{-1/2},length=200) 


cvr<-rep(0,length(lambdav))
for(i in 1:length(lambdav))
{
cvi=0
for(j in 1:5)
{
ind_pred<-indt[((j-1)*m+1):((j)*m)]
xpred=x[ind_pred,]
xest=x[-ind_pred,]
yest=y[-ind_pred]
ypred=y[ind_pred]
ttest=tvar[-ind_pred]
tpred=tvar[ind_pred]
result<-splinevcf(as.matrix(xest),yest,ttest,N,p,tpred,xpred,lambdav[i],0,weight,flag)
yest<-result[[3]]
cvi=cvi+mean((ypred-yest)^2)
print(mean((ypred-yest)^2))
}
cvr[i]=cvi
}
indx<-which.min(cvr)
lambda=lambdav[indx]
return(list(lambda,lambdav,cvr))
}







