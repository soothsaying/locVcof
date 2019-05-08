####################################################
#### Local varying coefficient model using the spline method in Xue, et al.
#### input: x is a nxd matrix
####        y: response variable
#####       tvar: a nx1 variable in the coefficient
#####       tpred: a npred x 1 vector, where the coefficients are evaluated
#####       xpred: a npred x d matrix, together with tpred, where yred are evaluated
#####       flag: 1: LASSO 2:SCAD 3:Adaptive LASSO 4: TLP 
####################################################
psvcf<-function(x,y,tvar,N,p,tpred,xpred)
{
n=length(y)
d=ncol(x)
temp=c(tvar,tpred)
###########################################################################
#### define the interior knots sequence: equally spaced knots (equal quantile knots can be
#### used as well
############################################################################
knots_n<-seq(min(tvar),max(tvar),length=N+2)  
knots_n<-knots_n[-c(1,(N+2))]
bsbasist<-bs(tvar,knots=knots_n,degree=p,intercept=TRUE,Boundary.knots=range(temp))
bsbasistpred<-bs(tpred,knots=knots_n,degree=p,intercept=TRUE,Boundary.knots=range(temp))

Jn=N+p+1
xdgn<-matrix(0,n,d*Jn)
for(j in 1:d)
{
temp<-((j-1)*Jn+1):(j*Jn)
xdgn[,temp]=bsbasist*matrix(rep(x[,j],Jn),ncol=Jn)
}

########Initial estimator the coefficients###
betahat<-solve(t(xdgn)%*%xdgn)%*%t(xdgn)%*%y
yhat<-xdgn%*%betahat
alphahat<-matrix(0,length(tpred),d)
for(j in 1:d)
{
indl=((j-1)*Jn+1):(j*Jn)
alphahat[,j]<-bsbasistpred%*%betahat[indl]
}
return(list(betahat,alphahat,yhat))
}