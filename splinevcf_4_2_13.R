####### Example with univariate variable which constrains certain values zero########
##
###N=number of interior knots
library(splines)
#library(glmnet)


############################################################
### Defines the first order derivative of LASSO,SCAD,TLP function
#############################################################
SCADder<-function(x,lambda,a)
{
absx=abs(x)
return(lambda*((absx<=lambda)+(absx>lambda)*(absx<=a*lambda)*((a*lambda-absx)/((a-1)*lambda))))
}

LASSOder<-function(x,lambda,weight)
{
absx=abs(x)
return(lambda*weight)
}

TLPder<-function(x,lambda,tau)
{
absx=abs(x)
return(lambda/tau*(absx/tau<1))
}

####################################################
#### Local varying coefficient model using the spline method in Xue, et al.
#### input: x is a nxd matrix
####        y: response variable
#####       tvar: a nx1 variable in the coefficient
#####       tpred: a npred x 1 vector, where the coefficients are evaluated
#####       xpred: a npred x d matrix, together with tpred, where yred are evaluated
#####       flag: 1: LASSO 2:SCAD 3:Adaptive LASSO 4: TLP 
####################################################
splinevcf<-function(x,y,tvar,N,p,tpred,xpred,lambda,tau,weight,flag)
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
betaini<-solve(t(xdgn)%*%xdgn)%*%t(xdgn)%*%y

Jn=N+p+1
betaold<-betaini
betahat<-betaold
epsilon=10^(-6)
diff<-1
step=0

while((diff>=epsilon)&(step<=100))
{

betaoldm<-matrix(betaold,ncol=d) #####put betaold into a matrix Jn x d
betaolds<-matrix(betaini,ncol=d) #####put betaold into a matrix Jn x d
cvec<-matrix(0,N+1,d) ####stores the "c_lj" constants defined in paper
for(i in 1:d)
{
for(j in 1:(N+1))
{
ind<-(j):(j+p)
jnorm=sqrt(sum((betaoldm[ind,i])^2))
jnorms=sqrt(sum((betaolds[ind,i])^2))
if(jnorm^2<=epsilon) betahat[(Jn*(i-1)+j):(Jn*(i-1)+j+p)]=0   ### Shrink betahat to zero if jnorm is small
else{ 
  if(flag==2){  ##SCAD
    a=3.7
    cvec[j,i]<-SCADder(jnorm,lambda,a)/jnorm
  }
  else if(flag==4)  
#    cvec[j,i]<-TLPder(jnorms,lambda,tau)/jnorm  ##TLP
    cvec[j,i]<-TLPder(jnorm,lambda,tau)/jnorm  ##TLP
  else             
    cvec[j,i]<-LASSOder(jnorm,lambda,weight[j,i])/jnorm  ##LASSO & Adaptive LASSO
}
#print(jnorm)
}
}

#####Update the nonzero components only in the following
###########################################################
Pn=length(betahat)
inds<-(1:Pn)[betahat!=0] ####Non zero components index

if(sum(betahat!=0)>=1)
{
xdgnsg<-as.matrix(xdgn[,inds])
sigma<-rep(0,length(inds))
cvec<-rbind(matrix(0,p,d),cvec,matrix(0,p,d))

indj=1
for(j in inds)
{
j1=j%%Jn
if(j1==0) j1=Jn
j2=ceiling(j/Jn)
sigma[indj]=sum(cvec[(j1):(j1+p),j2])
indj=indj+1
}

#print(inds)
if(length(inds)>1) SIGMA<-diag(sigma)
if(length(inds)<=1) {SIGMA=sigma}
betahat[inds]<-solve(t(xdgnsg)%*%xdgnsg+n*SIGMA)%*%t(xdgnsg)%*%y
diff=mean((betahat-betaold)^2)
step=step+1
betaold<-betahat
}
if(sum(betahat!=0)==0) diff=0
#print(diff)
#print(betahat)
}##end of while

alphahat<-matrix(0,length(tpred),d)
for(l in 1:d)
{
indl=((l-1)*Jn+1):(l*Jn)
alphahat[,l]<-bsbasistpred%*%betahat[indl]
}

ypred<-apply(alphahat*xpred,1,sum)

return(list(betahat,alphahat,ypred))
} ##end of function


############BIC tunning

lambdasel_BIC<-function(x,y,tvar,N,p,weight,flag)
{
#lambdav<-seq(0.0001,0.2,length=200)
n<-length(y)
m<-ceiling(n/5)
indt<-sample(n)

lambdav<-seq(0.01*(n/5)^{-1/5},1.5*(n/5)^{-1/5},length=200) #LASSO


###TLP
if(flag==4){
tauv<-seq(0.001,2,length=50)
BIC<-matrix(0,length(lambdav),length(tauv))
for(i in 1:length(lambdav))
{
for(j in 1:length(tauv))
{
result<-splinevcf(x,y,tvar,N,p,tvar,x,lambdav[i],tauv[j],0,flag)
betae<-result[[1]]
parp=sum(betae!=0)
yhat<-result[[3]]
BIC[i,j]=log(mean((y-yhat)^2))+log(n)*parp/n;
}
}
indx<-which(BIC==min(BIC),arr.ind=T)
lambda=lambdav[indx[1,1]]
tau=tauv[indx[1,2]]
return(list(lambda,tau,lambdav,tauv,BIC))
}

else{  ##SCAD LASSO & Adaptive LASSO

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


