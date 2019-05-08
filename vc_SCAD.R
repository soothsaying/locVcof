source("gSCAD.R")
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

vc_SCAD<-function(x,tvar,y,xe,te,N,p,lambda,a)
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
result_SCAD<-groupSCAD(dmatrix,y,NN,lambda,a)
####Prepare for outputs
spline_te<-bspline(te,knot,p)
alphahat_SCAD<-matrix(0,m,d)
sel_SCAD<-rep(0,d)

for(i in 1:d)
{
ind<-((i-1)*(N+p+1)+1):(i*(N+p+1))
alphahat_SCAD[,i]<-spline_te%*%result_SCAD[ind]
sel_SCAD[i]=1*(sum(abs(result_SCAD[ind]))!=0)
}
yhat_SCAD<-apply(alphahat_SCAD*xe,1,sum)

return(list(alphahat_SCAD,sel_SCAD,yhat_SCAD,result_SCAD))
}