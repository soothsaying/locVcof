source("gLasso.R")
####High dimensional varying coefficient model####
#######################################################
bspline<-function(x,knots,p)
{
N<-length(knots)
n<-length(x)
temp<-matrix(rep(x,N),ncol=N)-matrix(rep(knots,n*rep(1,N)),ncol=N) 
xmatrix<-(temp+abs(temp))/2
basis<-matrix(1,n,(1+N+p))
for(i in 1:p)
{
basis[,(i+1)]=x^i
}
basis[,((p+2):(N+p+1))]=xmatrix^p
return(basis)
}


##N is the number of interior knots

vc_LASSO<-function(x,tvar,y,xe,te,N,p,lambda,weight)
{
n<-length(y)
d<-dim(x)[[2]]
m<-length(te)

knot<-seq(min(tvar),max(tvar),length=(N+2))
knot=knot[-c(1,N+2)]
spline_t<-bspline(tvar,knot,p)
dmatrix<-matrix(0,n,(N+p+1)*d)
for(i in 1:d)
{
ind=((i-1)*(N+p+1)+1):(i*(N+p+1))
dmatrix[,ind]=spline_t*x[,i]
}
NN=N+p+1
result_LASSO<-groupLASSO(dmatrix,y,NN,lambda,weight)
####Prepare for outputs
spline_te<-bspline(te,knot,p)
alphahat_LASSO<-matrix(0,m,d)
sel_LASSO<-rep(0,d)

for(i in 1:d)
{
ind<-((i-1)*(N+p+1)+1):(i*(N+p+1))
alphahat_LASSO[,i]<-spline_te%*%result_LASSO[ind]
sel_LASSO[i]=1*(sum(abs(result_LASSO[ind]))!=0)
}
yhat_LASSO<-apply(alphahat_LASSO*xe,1,sum)

return(list(alphahat_LASSO,sel_LASSO,yhat_LASSO,result_LASSO))
}