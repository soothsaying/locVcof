library("mvtnorm")
library("locVcof")
#######################################
###Test Data
########################################
n=500  ####number of subject
m=5    #### number of observation per subject
p=2    ### Order of polynomial spline
nt=n*m
N=2*ceiling(nt^{1/(2*p+3)}) # number of spline knots
rho=0.5   ###Moderate correlation 
d=5 ### number of independent variables

#case 1 
testf1<-function(x) return(rep(0,length(x)))
#case 2 
testf2<-function(x) sin(2*pi*x)
#case 3 
testf3<-function(x) 5*(((x-0.5)^2-0.025)*(x<0.342)+(-(x-0.5)^2+0.025)*(x>0.658))
#case 4 
testf4<-function(x) 3*((-x+0.3)*(x<0.3)+(x-0.7)*(x>0.7))
#case 5 
testf5<-function(x) (-45/2*(x-0.5)^2+0.9)*(x>0.3&x<0.7)



set.seed(1234)


meanep<-rep(0,m)
sigmaep<-matrix(0.5*rho,m,m)+diag(0.5-0.5*rho,m,m)

# generate AR(1) data
f_ar1<-function(epsi, a){
  epsiar = epsi
  for (i in 2:dim(epsi)[1]){
    epsiar[i,] = a*epsiar[i-1,]+epsiar[i,]
  }
  epsiar 
}


x<-matrix(rnorm(nt*d,0,1),ncol=d)
tvar<-runif(nt,0,1)
tpred<-seq(0,1,length=5000)
xpred<-matrix(rnorm(d*5000,0,1),5000,d)
epsi<-rmvnorm(n,mean=meanep,sigma=sigmaep)
y<-testf1(tvar)*x[,1]+testf2(tvar)*x[,2]+testf3(tvar)*x[,3]+testf4(tvar)*x[,4]+testf5(tvar)*x[,5]+c(t(f_ar1(epsi, 0.5)))

alpha_TRUE=cbind(testf1(tpred),testf2(tpred),testf3(tpred),testf4(tpred),testf5(tpred))
##################Full estimation ###############################
result_full<-psvcf(x,y,tvar,N,p,tpred,xpred)
alphahat<-result_full[[2]]



###SCAD
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,0,flag=2)
lambda_SCAD<-lambdaresult[[1]]
result_SCAD<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_SCAD,0,0,flag=2)
alphahat_SCAD<-result_SCAD[[2]]
betahat_SCAD<-matrix(result_SCAD[[1]],nrow=N+p+1,ncol=d)




###TLP
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,0,flag=4)
lambda_TLP<-lambdaresult[[1]]
tau_TLP<-lambdaresult[[2]]
result_TLP<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_TLP,tau_TLP,0,flag=4)
alphahat_TLP<-result_TLP[[2]]
betahat_TLP<-matrix(result_TLP[[1]],nrow=N+p+1,ncol=d)


###LASSO
weight<-matrix(1,nrow=N+1,ncol=d)
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,weight,flag=1)
lambda<-lambdaresult[[1]]
result_LASSO<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda,0,weight,flag=1)
alphahat_LASSO<-result_LASSO[[2]]
betahat_LASSO<-matrix(result_LASSO[[1]],nrow=N+p+1,ncol=d)



###Adaptive LASSO
result_LSE<-splinevcf(x,y,tvar,N,p,tpred,xpred,0,0,weight,flag=1)
alphahat_LSE<-result_LSE[[2]]
betahat_LSE<-matrix(result_LSE[[1]],nrow=N+p+1,ncol=d)
weight<-matrix(10^8,nrow=N+1,ncol=d)

for (i in 1:d){
   for(j in 1:(N+1)){
       temp=sum((betahat_LSE[j:(j+p),i])^2)
       if(temp!=0)
          weight[j,i]=1/temp
   }
}
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,weight,flag=3)
lambda<-lambdaresult[[1]]
result_ADLASSO<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda,0,weight,flag=3)
alphahat_ADLASSO<-result_ADLASSO[[2]]
