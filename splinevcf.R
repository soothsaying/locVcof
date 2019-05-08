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
cvec<-matrix(0,N+1,d) ####stores the "c_lj" constants defined in paper
for(i in 1:d)
{
for(j in 1:(N+1))
{
ind<-(j):(j+p)
jnorm=sqrt(sum((betaoldm[ind,i])^2))
if(jnorm^2<=epsilon) betahat[(Jn*(i-1)+j):(Jn*(i-1)+j+p)]=0   ### Shrink betahat to zero if jnorm is small
else{ 
  if(flag==2){  ##SCAD
    a=3.7
    cvec[j,i]<-SCADder(jnorm,lambda,a)/jnorm
  }
  else if(flag==4)  
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


############cross validation function
lambdasel<-function(x,y,tvar,N,p,weight,flag)
{
lambdav<-seq(0.0001,0.2,length=200)
n<-length(y)
m<-ceiling(n/5)
indt<-sample(n)



###TLP
if(flag==4){
tauv<-seq(0.0001,0.5,length=50)
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

else{  ##SCAD LASSO & Adaptive LASSO

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


}



#######################################
###Test Data
########################################
n=400
N=5
p=2



#case 1 beta=0 N=3
testf1<-function(x) return(rep(0,length(x)))
#case 2 N=5 p=2
testf2<-function(x) sin(2*pi*x)
#case 3 N=5 p=2
testf3<-function(x) 5*(((x-0.5)^2-0.025)*(x<0.342)+(-(x-0.5)^2+0.025)*(x>0.658))
#case 4 N=5 p=2
testf4<-function(x) 3*((-x+0.3)*(x<0.3)+(x-0.7)*(x>0.7))
#case 5 N=9 p=2
testf5<-function(x) (-45/2*(x-0.5)^2+0.9)*(x>0.3&x<0.7)
#case 6 N=10 p=2
testf6<-function(x){
   len<-length(x)
   fx<-rep(0, len)
   for(i in 1:len)
   { 
      if(x[i]>0.25 && x[i]<0.75) 
           fx[i]=sin(2*pi*x[i])
      else
           fx[i]=0
   }
   return(fx)
}
#case 7 N=15 p=2
testf7<-function(x) -sin(5*pi*x)*(x<0.4&x>0.2)+sin(5*pi*x)*(x>0.6&x<0.8)


d=1
x<-matrix(rnorm(n*d,0,1),ncol=d)
tvar<-runif(n,0,1)
tpred<-seq(0,1,length=100)
xpred<-as.matrix(x[1:length(tpred),])
eps<-rnorm(100,0,1)

FLAG=3 # represents num of case
y<-testf3(tvar)*x[,1]+0.5*rnorm(n,0,1)
ypred<-testf3(tpred)*xpred[,1]+0.5*eps
alpha_TRUE=testf3(tpred)
#alpha_TRUE=cbind(testf2(tpred),matrix(0,nrow=100,ncol=d-4))


###LASSO
d=dim(x)[2]
weight<-matrix(1,nrow=N+1,ncol=d)
lambdaresult<-lambdasel(x,y,tvar,N,p,weight,flag=1)
lambda_LASSO<-lambdaresult[[1]]
result_LASSO<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_LASSO,0,weight,flag=1)
alphahat_LASSO<-result_LASSO[[2]]
betahat_LASSO<-matrix(result_LASSO[[1]],nrow=N+p+1,ncol=d)

###SCAD
lambdaresult<-lambdasel(x,y,tvar,N,p,0,flag=2)
lambda_SCAD<-lambdaresult[[1]]
result_SCAD<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_SCAD,0,0,flag=2)
alphahat_SCAD<-result_SCAD[[2]]
betahat_SCAD<-matrix(result_SCAD[[1]],nrow=N+p+1,ncol=d)

###Adaptive LASSO
weight<-matrix(10^8,nrow=N+1,ncol=d)

for (i in 1:d){
   for(j in 1:(N+1)){
       temp=sum((betahat_SCAD[j:(j+p),i])^2)
       if(temp!=0)
          weight[j,i]=1/temp
   }
}
lambdaresult<-lambdasel(x,y,tvar,N,p,weight,flag=3)
lambda_AdLASSO<-lambdaresult[[1]]
result_AdLASSO<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_AdLASSO,0,weight,flag=3)
alphahat_AdLASSO<-result_AdLASSO[[2]]

###TLP
lambdaresult<-lambdasel(x,y,tvar,N,p,0,flag=4)
lambda_TLP<-lambdaresult[[1]]
tau_TLP<-lambdaresult[[2]]
result_TLP<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_TLP,tau_TLP,0,flag=4)
alphahat_TLP<-result_TLP[[2]]


postscript(file="case6a.ps",height=2,width=6)
layout(matrix(c(1,0,2),1,3),c(3,0,3),2)
par(mar=c(5,5,2,2)+.1,mex=.6)
plot(tpred,alpha_TRUE,col=1,type='l',lwd=2,ylim=c(-1,1),xlab="T",ylab="Alpha")
dev.off()

postscript(file="case6b.ps",height=2,width=6)
layout(matrix(c(1,0,2),1,3),c(3,0,3),2)
par(mar=c(5,5,2,2)+.1,mex=.6)
plot(tpred,alpha_TRUE,col=1,type='l',ylim=c(-1,1),xlab="T",ylab="Alpha")
lines(tpred,alphahat_SCAD[,1],lty=2,col=2,lwd=2)
lines(tpred,alphahat_LASSO[,1],lty=3,col=3,lwd=2)
lines(tpred,alphahat_AdLASSO[,1],lty=4,col=4,lwd=2)
lines(tpred,alphahat_TLP[,1],lty=5,col=5,lwd=2)
legtxt<-c(expression("True"),expression("SCAD"),expression("LASSO"),expression("AdLASSO"),expression("TLP"))
legend(0.6,1.05,col=1:5,lty=1:5,legend=legtxt,lwd=c(2,2,2,2,2),cex=0.8)
dev.off()

