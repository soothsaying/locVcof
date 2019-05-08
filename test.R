source("splinevcf_4_2_13_two_int.R")
source("psvcf.R")
#######################################
###Test Data
########################################
n=500
p=1
N=2*floor(n^{1/(2*p+3)})


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


set.seed(1234)
d=5
nrep=100

mse_SCAD<-matrix(0,nrep,d)
max_SCAD<-matrix(0,nrep,d)
Null_SCAD<-matrix(0,nrep,d)
Sig_SCAD<-matrix(0,nrep,d)

mse_TLP<-matrix(0,nrep,d)
max_TLP<-matrix(0,nrep,d)
Null_TLP<-matrix(0,nrep,d)
Sig_TLP<-matrix(0,nrep,d)

mse_LASSO<-matrix(0,nrep,d)
max_LASSO<-matrix(0,nrep,d)
Null_LASSO<-matrix(0,nrep,d)
Sig_LASSO<-matrix(0,nrep,d)

mse_ADLASSO<-matrix(0,nrep,d)
max_ADLASSO<-matrix(0,nrep,d)
Null_ADLASSO<-matrix(0,nrep,d)
Sig_ADLASSO<-matrix(0,nrep,d)


int_SCAD<-matrix(0,nrep,2)
int_TLP<-matrix(0,nrep,2)
int_LASSO<-matrix(0,nrep,2)
int_ADLASSO<-matrix(0,nrep,2)


for(r in 1:2)
{
x<-matrix(rnorm(n*d,0,1),ncol=d)
tvar<-runif(n,0,1)
tpred<-seq(0,1,length=5000)
xpred<-matrix(rnorm(d*5000,0,1),5000,d)
y<-testf1(tvar)*x[,1]+testf2(tvar)*x[,2]+testf3(tvar)*x[,3]+testf4(tvar)*x[,4]+testf5(tvar)*x[,5]+0.5*rnorm(n,0,1)

alpha_TRUE=cbind(testf1(tpred),testf2(tpred),testf3(tpred),testf4(tpred),testf5(tpred))
##################Full estimation ###############################
result_full<-psvcf(x,y,tvar,N,p,tpred,xpred)
alphahat<-result_full[[2]]
mse_FULL[r,]<-apply((alphahat-alpha_TRUE)^2,2,mean)
max_FULL[r,]<-apply(abs(alphahat-alpha_TRUE),2,max)


###SCAD
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,0,flag=2)
lambda_SCAD<-lambdaresult[[1]]
result_SCAD<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_SCAD,0,0,flag=2)
alphahat_SCAD<-result_SCAD[[2]]
betahat_SCAD<-matrix(result_SCAD[[1]],nrow=N+p+1,ncol=d)
mse_SCAD[r,]=apply((alphahat_SCAD-alpha_TRUE)^2,2,mean)
max_SCAD[r,]<-apply(abs(alphahat_SCAD-alpha_TRUE),2,max)

ind1=(1:5000)
ind2=(1:5000)
ind3=(1:5000)[(tpred>=0.342)&(tpred<0.658)]
ind4=(1:5000)[(tpred>=0.3)&(tpred<=0.7)]
ind5=(1:5000)[(tpred<=0.3)|(tpred>=0.7)]

TR=alpha_TRUE==0
SCAD=alphahat_SCAD==0

Null_SCAD[r,1]=mean(SCAD[ind1,1]==TR[ind1,1])
Null_SCAD[r,3]=mean(SCAD[ind3,3]==TR[ind3,3])
Null_SCAD[r,4]=mean(SCAD[ind4,4]==TR[ind4,4])
Null_SCAD[r,5]=mean(SCAD[ind5,5]==TR[ind5,5])

Sig_SCAD[r,2]=mean(SCAD[ind2,2]==TR[ind2,2])
Sig_SCAD[r,3]=mean(SCAD[-ind3,3]==TR[-ind3,3])
Sig_SCAD[r,4]=mean(SCAD[-ind4,4]==TR[-ind4,4])
Sig_SCAD[r,5]=mean(SCAD[-ind5,5]==TR[-ind5,5])

int_SCAD[r,1]<-min(tpred[SCAD[,3]])
int_SCAD[r,2]<-max(tpred[SCAD[,3]])



###TLP
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,0,flag=4)
lambda_TLP<-lambdaresult[[1]]
tau_TLP<-lambdaresult[[2]]
result_TLP<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_TLP,tau_TLP,0,flag=4)
alphahat_TLP<-result_TLP[[2]]

mse_TLP[r,]=apply((alphahat_TLP-alpha_TRUE)^2,2,mean)
max_TLP[r,]<-apply(abs(alphahat_TLP-alpha_TRUE),2,max)


TLP=alphahat_TLP==0
Null_TLP[r,1]=mean(TLP[ind1,1]==TR[ind1,1])
Null_TLP[r,3]=mean(TLP[ind3,3]==TR[ind3,3])
Null_TLP[r,4]=mean(TLP[ind4,4]==TR[ind4,4])
Null_TLP[r,5]=mean(TLP[ind5,5]==TR[ind5,5])

Sig_TLP[r,2]=mean(TLP[ind2,2]==TR[ind2,2])
Sig_TLP[r,3]=mean(TLP[-ind3,3]==TR[-ind3,3])
Sig_TLP[r,4]=mean(TLP[-ind4,4]==TR[-ind4,4])
Sig_TLP[r,5]=mean(TLP[-ind5,5]==TR[-ind5,5])

int_TLP[r,1]<-min(tpred[TLP[,3]])
int_TLP[r,2]<-max(tpred[TLP[,3]])



###LASSO
weight<-matrix(1,nrow=N+1,ncol=d)
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,weight,flag=1)
lambda<-lambdaresult[[1]]
result_LASSO<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda,0,weight,flag=1)
alphahat_LASSO<-result_LASSO[[2]]
betahat_LASSO<-matrix(result_LASSO[[1]],nrow=N+p+1,ncol=d)
alphahat_LASSO[abs(alphahat_LASSO)<10^-3]=0

mse_LASSO[r,]=apply((alphahat_LASSO-alpha_TRUE)^2,2,mean)
max_LASSO[r,]<-apply(abs(alphahat_LASSO-alpha_TRUE),2,max)


LASSO=alphahat_LASSO==0
Null_LASSO[r,1]=mean(LASSO[ind1,1]==TR[ind1,1])
Null_LASSO[r,3]=mean(LASSO[ind3,3]==TR[ind3,3])
Null_LASSO[r,4]=mean(LASSO[ind4,4]==TR[ind4,4])
Null_LASSO[r,5]=mean(LASSO[ind5,5]==TR[ind5,5])

Sig_LASSO[r,2]=mean(LASSO[ind2,2]==TR[ind2,2])
Sig_LASSO[r,3]=mean(LASSO[-ind3,3]==TR[-ind3,3])
Sig_LASSO[r,4]=mean(LASSO[-ind4,4]==TR[-ind4,4])
Sig_LASSO[r,5]=mean(LASSO[-ind5,5]==TR[-ind5,5])

int_LASSO[r,1]<-min(tpred[LASSO[,3]])
int_LASSO[r,2]<-max(tpred[LASSO[,3]])

###Adaptive LASSO
weight<-matrix(10^8,nrow=N+1,ncol=d)

for (i in 1:d){
   for(j in 1:(N+1)){
       temp=sum((betahat_LASSO[j:(j+p),i])^2)
       if(temp!=0)
          weight[j,i]=1/temp
   }
}
lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,weight,flag=3)
lambda<-lambdaresult[[1]]
result_ADLASSO<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda,0,weight,flag=3)
alphahat_ADLASSO<-result_ADLASSO[[2]]
#alphahat_AdLASSO[abs(alphahat_ADLASSO)<10^-3]=0

mse_ADLASSO[r,]=apply((alphahat_ADLASSO-alpha_TRUE)^2,2,mean)
max_ADLASSO[r,]<-apply(abs(alphahat_ADLASSO-alpha_TRUE),2,max)


ADLASSO=alphahat_ADLASSO==0
Null_ADLASSO[r,1]=mean(ADLASSO[ind1,1]==TR[ind1,1])
Null_ADLASSO[r,3]=mean(ADLASSO[ind3,3]==TR[ind3,3])
Null_ADLASSO[r,4]=mean(ADLASSO[ind4,4]==TR[ind4,4])
Null_ADLASSO[r,5]=mean(ADLASSO[ind5,5]==TR[ind5,5])

Sig_ADLASSO[r,2]=mean(ADLASSO[ind2,2]==TR[ind2,2])
Sig_ADLASSO[r,3]=mean(ADLASSO[-ind3,3]==TR[-ind3,3])
Sig_ADLASSO[r,4]=mean(ADLASSO[-ind4,4]==TR[-ind4,4])
Sig_ADLASSO[r,5]=mean(ADLASSO[-ind5,5]==TR[-ind5,5])

int_ADLASSO[r,1]<-min(tpred[ADLASSO[,3]])
int_ADLASSO[r,2]<-max(tpred[ADLASSO[,3]])


print(r)
}


mean(apply(Null_SCAD[,-2],2,mean))
mean(apply(Null_TLP[,-2],2,mean))
mean(apply(Null_LASSO[,-2],2,mean))
mean(apply(Null_ADLASSO[,-2],2,mean))

mean(apply(Sig_SCAD[,-1],2,mean))
mean(apply(Sig_TLP[,-1],2,mean))
mean(apply(Sig_LASSO[,-1],2,mean))
mean(apply(Sig_ADLASSO[,-1],2,mean))


mean(apply(mse_TLP,2,mean))
mean(apply(mse_SCAD,2,mean))
mean(apply(mse_LASSO,2,mean))
mean(apply(mse_ADLASSO,2,mean))


mean(apply(max_TLP,2,mean))
mean(apply(max_SCAD,2,mean))
mean(apply(max_LASSO,2,mean))
mean(apply(max_ADLASSO,2,mean))

#int_TLP[c(5,67),1]=rep(0,1)
#int_TLP[c(5,67),2]=rep(1,1)

apply(int_TLP,2,mean)
apply(int_SCAD,2,mean)
apply(int_LASSO,2,mean)
apply(int_ADLASSO,2,mean)


par(mfrow=c(2,3))
plot(tpred,testf1(tpred),type='l')
lines(tpred,alphahat_SCAD[,1],col=2)
lines(tpred,alphahat_TLP[,1],col=3)

plot(tpred,testf2(tpred),type='l',ylim=c(-1.5,1.5))
lines(tpred,alphahat_SCAD[,2],col=2)
lines(tpred,alphahat_TLP[,2],col=3)


plot(tpred,testf3(tpred),type='l')
lines(tpred,alphahat_SCAD[,3],col=2)
lines(tpred,alphahat_TLP[,3],col=3)

plot(tpred,testf4(tpred),type='l',ylim=c(-0.5,1))
lines(tpred,alphahat_SCAD[,4],col=2)
lines(tpred,alphahat_TLP[,4],col=3)

plot(tpred,testf5(tpred),type='l')
lines(tpred,alphahat_SCAD[,5],col=2)
lines(tpred,alphahat_TLP[,5],col=3)



