setwd("C:/Users/wdsxsx/Desktop/R_code")
library("splines")
library("mvtnorm")
source("splinevcf_4_2_13_two_int.R")
source("psvcf.R")
library(randomcoloR)
#######################################
###Test Data
########################################
n=100  ####number of subject
m=5    #### number of observation per subject
p=1    ### Order of polynomial spline
nt=n*m
N=2*ceiling(nt^{1/(2*p+3)})
rho=0.5   ###Moderate correlation 


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


x<-matrix(rnorm(nt*d,0,1),ncol=d)
tvar<-runif(nt,0,1)
tpred<-seq(0,1,length=5000)
xpred<-matrix(rnorm(d*5000,0,1),5000,d)

meanep<-rep(0,m)
sigmaep<-matrix(0.5*rho,m,m)+diag(0.5-0.5*rho,m,m)
epsi<-rmvnorm(n,mean=meanep,sigma=sigmaep)
y<-testf1(tvar)*x[,1]+testf2(tvar)*x[,2]+testf3(tvar)*x[,3]+testf4(tvar)*x[,4]+testf5(tvar)*x[,5]+c(t(epsi))

alpha_TRUE=cbind(testf1(tpred),testf2(tpred),testf3(tpred),testf4(tpred),testf5(tpred))
TR=alpha_TRUE==0

lambdaresult<-lambdasel_BIC(x,y,tvar,N,p,0,flag=4)
lambda_TLP<-lambdaresult[[1]]
tau_TLP<-lambdaresult[[2]]

lambdav<-seq(0.001*(n/5)^{-1/5}, 0.2*(n/5)^{-1/5},length=200)
tauv<-seq(0.001,2,length=50)

beta_lam_testf1<-matrix(0,length(lambdav), N+p+1)
beta_lam_testf2<-matrix(0,length(lambdav), N+p+1)
beta_lam_testf3<-matrix(0,length(lambdav), N+p+1)
beta_lam_testf4<-matrix(0,length(lambdav), N+p+1)
beta_lam_testf5<-matrix(0,length(lambdav), N+p+1)

beta_tau_testf1<-matrix(0,length(tauv), N+p+1)
beta_tau_testf2<-matrix(0,length(tauv), N+p+1)
beta_tau_testf3<-matrix(0,length(tauv), N+p+1)
beta_tau_testf4<-matrix(0,length(tauv), N+p+1)
beta_tau_testf5<-matrix(0,length(tauv), N+p+1)

alpha_lam_testf1<-matrix(0,length(lambdav), N+1)
alpha_lam_testf2<-matrix(0,length(lambdav), N+1)
alpha_lam_testf3<-matrix(0,length(lambdav), N+1)
alpha_lam_testf4<-matrix(0,length(lambdav), N+1)
alpha_lam_testf5<-matrix(0,length(lambdav), N+1)

alpha_tau_testf1<-matrix(0,length(tauv), N+1)
alpha_tau_testf2<-matrix(0,length(tauv), N+1)
alpha_tau_testf3<-matrix(0,length(tauv), N+1)
alpha_tau_testf4<-matrix(0,length(tauv), N+1)
alpha_tau_testf5<-matrix(0,length(tauv), N+1)

Sig_TLP_lam<-matrix(0,length(lambdav), d)
Sig_TLP_tau<-matrix(0,length(tauv), d)

Null_TLP_lam<-matrix(0,length(lambdav), d)
Null_TLP_tau<-matrix(0,length(tauv), d)

ind1=(1:5000)
ind2=(1:5000)
ind3=(1:5000)[(tpred>=0.342)&(tpred<0.658)]
ind4=(1:5000)[(tpred>=0.3)&(tpred<=0.7)]
ind5=(1:5000)[(tpred<=0.3)|(tpred>=0.7)]


for(i in 1:length(lambdav)){
  result<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambdav[i],tau_TLP,0,flag=4)
  beta = matrix(result[[1]], N+p+1, d)
  beta_lam_testf1[i,] <- beta[,1]
  beta_lam_testf2[i,] <- beta[,2]
  beta_lam_testf3[i,] <- beta[,3]
  beta_lam_testf4[i,] <- beta[,4]  
  beta_lam_testf5[i,] <- beta[,5]
  for(j in 1:(N+1)){
    alpha_lam_testf1[i,j] <- sqrt(beta_lam_testf1[i, j]**2 + beta_lam_testf1[i, j+1]**2)
    alpha_lam_testf2[i,j] <- sqrt(beta_lam_testf2[i, j]**2 + beta_lam_testf2[i, j+1]**2)
    alpha_lam_testf3[i,j] <- sqrt(beta_lam_testf3[i, j]**2 + beta_lam_testf3[i, j+1]**2)
    alpha_lam_testf4[i,j] <- sqrt(beta_lam_testf4[i, j]**2 + beta_lam_testf4[i, j+1]**2)
    alpha_lam_testf5[i,j] <- sqrt(beta_lam_testf5[i, j]**2 + beta_lam_testf5[i, j+1]**2)
  }

  
  TLP=result[[2]]==0

  
  Sig_TLP_lam[i,2]=mean(TLP[ind2,2]==TR[ind2,2])
  Sig_TLP_lam[i,3]=mean(TLP[-ind3,3]==TR[-ind3,3])
  Sig_TLP_lam[i,4]=mean(TLP[-ind4,4]==TR[-ind4,4])
  Sig_TLP_lam[i,5]=mean(TLP[-ind5,5]==TR[-ind5,5])
  
  Null_TLP_lam[i,1]=mean(TLP[ind1,1]==TR[ind1,1])
  Null_TLP_lam[i,3]=mean(TLP[ind3,3]==TR[ind3,3])
  Null_TLP_lam[i,4]=mean(TLP[ind4,4]==TR[ind4,4])
  Null_TLP_lam[i,5]=mean(TLP[ind5,5]==TR[ind5,5])
}

for(i in 1:length(tauv)){
  result<-splinevcf(x,y,tvar,N,p,tpred,xpred,lambda_TLP,tauv[i],0,flag=4)
  beta = matrix(result[[1]], N+p+1, d)
  
  beta_tau_testf1[i,] <- beta[,1]
  beta_tau_testf2[i,] <- beta[,2]
  beta_tau_testf3[i,] <- beta[,3]
  beta_tau_testf4[i,] <- beta[,4]  
  beta_tau_testf5[i,] <- beta[,5]

  for(j in 1:(N+1)){
    alpha_tau_testf1[i,j] <- sqrt(beta_tau_testf1[i, j]**2 + beta_tau_testf1[i, j+1]**2)
    alpha_tau_testf2[i,j] <- sqrt(beta_tau_testf2[i, j]**2 + beta_tau_testf2[i, j+1]**2)
    alpha_tau_testf3[i,j] <- sqrt(beta_tau_testf3[i, j]**2 + beta_tau_testf3[i, j+1]**2)
    alpha_tau_testf4[i,j] <- sqrt(beta_tau_testf4[i, j]**2 + beta_tau_testf4[i, j+1]**2)
    alpha_tau_testf5[i,j] <- sqrt(beta_tau_testf5[i, j]**2 + beta_tau_testf5[i, j+1]**2)
  }
  
  
  
  TLP=result[[2]]==0
  

  Sig_TLP_tau[i,2]=mean(TLP[ind2,2]==TR[ind2,2])
  Sig_TLP_tau[i,3]=mean(TLP[-ind3,3]==TR[-ind3,3])
  Sig_TLP_tau[i,4]=mean(TLP[-ind4,4]==TR[-ind4,4])
  Sig_TLP_tau[i,5]=mean(TLP[-ind5,5]==TR[-ind5,5])
  
  
  Null_TLP_tau[i,1]=mean(TLP[ind1,1]==TR[ind1,1])
  Null_TLP_tau[i,3]=mean(TLP[ind3,3]==TR[ind3,3])
  Null_TLP_tau[i,4]=mean(TLP[ind4,4]==TR[ind4,4])
  Null_TLP_tau[i,5]=mean(TLP[ind5,5]==TR[ind5,5])
}


par(mfrow=c(2,2))
plot(lambdav,Sig_TLP_lam[,2],type='l',lwd=2, xlab=expression(lambda), ylab='SIG', main=expression(alpha[2](T)))
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.02, 1.2, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)
plot(lambdav,Sig_TLP_lam[,3],type='l',lwd=2, xlab=expression(lambda), ylab='SIG', main=expression(alpha[3](T)))
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.02, 0.95, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)
plot(lambdav,Sig_TLP_lam[,4],type='l',lwd=2, xlab=expression(lambda), ylab='SIG', main=expression(alpha[4](T)))
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.02, 0.95, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)
plot(lambdav,Sig_TLP_lam[,5],type='l',lwd=2, xlab=expression(lambda), ylab='SIG', main=expression(alpha[5](T)))
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.02, 0.98, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)


par(mfrow=c(2,2))
plot(tauv,Sig_TLP_tau[,2],type='l',lwd=2, xlab=expression(tau), ylab='SIG', main=expression(alpha[2](T)))
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(0.8, 1, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)
plot(tauv,Sig_TLP_tau[,3],type='l',lwd=2, xlab=expression(tau), ylab='SIG', main=expression(alpha[3](T)))
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(0.8, 1, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)
plot(tauv,Sig_TLP_tau[,4],type='l',lwd=2, xlab=expression(tau), ylab='SIG', main=expression(alpha[4](T)))
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(0.8, 1, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)
plot(tauv,Sig_TLP_tau[,5],type='l',lwd=2, xlab=expression(tau), ylab='SIG', main=expression(alpha[5](T)))
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(0.8, 1, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.7)


colors <- distinctColorPalette(N+1)

#### lambda plot
leg = c(expression(gamma[1][1]), expression(gamma[1][2]), expression(gamma[1][3]), expression(gamma[1][4]), expression(gamma[1][5]), 
        expression(gamma[1][6]), expression(gamma[1][7]), expression(gamma[1][8]), expression(gamma[1][9]))
par(mfrow=c(1,1))
plot(lambdav,alpha_lam_testf1[,1],type='l',lwd=2, col=colors[1], xlab=expression(lambda), ylab='L2 norm', 
     ylim=c(min(alpha_lam_testf1)-0.05, max(alpha_lam_testf1)+0.05), main=expression(alpha[1](T)))
for(i in 2:(N+1)){
  lines(lambdav,alpha_lam_testf1[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.1, max(alpha_lam_testf1)+0.05, legend=leg, lwd=2, col=colors)
legend(0.09, max(alpha_lam_testf1)-0.08, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)



leg = c(expression(gamma[2][1]), expression(gamma[2][2]), expression(gamma[2][3]), expression(gamma[2][4]), expression(gamma[2][5]), 
        expression(gamma[2][6]), expression(gamma[2][7]), expression(gamma[2][8]), expression(gamma[2][9]))
par(mfrow=c(1,1))
plot(lambdav,alpha_lam_testf2[,1],type='l',lwd=2, col=colors[1], xlab=expression(lambda), ylab='L2 norm', 
     ylim=c(min(alpha_lam_testf2)-0.05, max(alpha_lam_testf2)+0.05), main=expression(alpha[2](T)))
for(i in 2:(N+1)){
  lines(lambdav,alpha_lam_testf2[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.1, max(alpha_lam_testf2)+0.05, legend=leg, lwd=2, col=colors)
legend(0.09, max(alpha_lam_testf2)-0.3, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)



leg = c(expression(gamma[3][1]), expression(gamma[3][2]), expression(gamma[3][3]), expression(gamma[3][4]), expression(gamma[3][5]), 
        expression(gamma[3][6]), expression(gamma[3][7]), expression(gamma[3][8]), expression(gamma[3][9]))
par(mfrow=c(1,1))
plot(lambdav,alpha_lam_testf3[,1],type='l',lwd=2, col=colors[1], xlab=expression(lambda), ylab='L2 norm', 
     ylim=c(min(alpha_lam_testf3)-0.05, max(alpha_lam_testf3)+0.05), main=expression(alpha[3](T)))
for(i in 2:(N+1)){
  lines(lambdav,alpha_lam_testf3[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.1, max(alpha_lam_testf3)+0.05, legend=leg, lwd=2, col=colors)
legend(0.09, max(alpha_lam_testf3)-0.5, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)



leg = c(expression(gamma[4][1]), expression(gamma[4][2]), expression(gamma[4][3]), expression(gamma[4][4]), expression(gamma[4][5]), 
        expression(gamma[4][6]), expression(gamma[4][7]), expression(gamma[4][8]), expression(gamma[4][9]))
par(mfrow=c(1,1))
plot(lambdav,alpha_lam_testf4[,1],type='l',lwd=2, col=colors[1], xlab=expression(lambda), ylab='L2 norm', 
     ylim=c(min(alpha_lam_testf4)-0.05, max(alpha_lam_testf4)+0.05), main=expression(alpha[4](T)))
for(i in 2:(N+1)){
  lines(lambdav,alpha_lam_testf4[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.1, max(alpha_lam_testf4)+0.05, legend=leg, lwd=2, col=colors)
legend(0.09, max(alpha_lam_testf4)-0.4, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)




leg = c(expression(gamma[5][1]), expression(gamma[5][2]), expression(gamma[5][3]), expression(gamma[5][4]), expression(gamma[5][5]), 
        expression(gamma[5][6]), expression(gamma[5][7]), expression(gamma[5][8]), expression(gamma[5][9]))
par(mfrow=c(1,1))
plot(lambdav,alpha_lam_testf5[,1],type='l',lwd=2, col=colors[1], xlab=expression(lambda), ylab='L2 norm', 
     ylim=c(min(alpha_lam_testf5)-0.05, max(alpha_lam_testf5)+0.05), main=expression(alpha[5](T)))
for(i in 2:(N+1)){
  lines(lambdav,alpha_lam_testf5[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=lambda_TLP, col = "red", lty = 2, lwd=2)
legend(0.1, max(alpha_lam_testf5)+0.05, legend=leg, lwd=2, col=colors)
legend(0.09, max(alpha_lam_testf5)-0.5, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)


# tau plot

leg = c(expression(gamma[1][1]), expression(gamma[1][2]), expression(gamma[1][3]), expression(gamma[1][4]), expression(gamma[1][5]), 
        expression(gamma[1][6]), expression(gamma[1][7]), expression(gamma[1][8]), expression(gamma[1][9]))
par(mfrow=c(1,1))
plot(tauv,alpha_tau_testf1[,1],type='l',lwd=2, col=colors[1], xlab=expression(tau), ylab='L2 norm', 
     ylim=c(min(alpha_tau_testf1)-0.05, max(alpha_tau_testf1)+0.05), main=expression(alpha[1](T)))
for(i in 2:(N+1)){
  lines(tauv,alpha_tau_testf1[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(1.8, max(alpha_tau_testf1)+0.05, legend=leg, lwd=2, col=colors)
legend(1.6, max(alpha_tau_testf1)-0.1, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)


leg = c(expression(gamma[2][1]), expression(gamma[2][2]), expression(gamma[2][3]), expression(gamma[2][4]), expression(gamma[2][5]), 
        expression(gamma[2][6]), expression(gamma[2][7]), expression(gamma[2][8]), expression(gamma[2][9]))
par(mfrow=c(1,1))
plot(tauv,alpha_tau_testf2[,1],type='l',lwd=2, col=colors[1], xlab=expression(tau), ylab='L2 norm', 
     ylim=c(min(alpha_tau_testf2)-0.05, max(alpha_tau_testf2)+0.05), main=expression(alpha[2](T)))
for(i in 2:(N+1)){
  lines(tauv,alpha_tau_testf2[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(1.8, max(alpha_tau_testf2)+0.05, legend=leg, lwd=2, col=colors)
legend(1.6, max(alpha_tau_testf2)-0.6, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)



leg = c(expression(gamma[3][1]), expression(gamma[3][2]), expression(gamma[3][3]), expression(gamma[3][4]), expression(gamma[3][5]), 
        expression(gamma[3][6]), expression(gamma[3][7]), expression(gamma[3][8]), expression(gamma[3][9]))
par(mfrow=c(1,1))
plot(tauv,alpha_tau_testf3[,1],type='l',lwd=2, col=colors[1], xlab=expression(tau), ylab='L2 norm', 
     ylim=c(min(alpha_tau_testf3)-0.05, max(alpha_tau_testf3)+0.05), main=expression(alpha[3](T)))
for(i in 2:(N+1)){
  lines(tauv,alpha_tau_testf3[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(1.8, max(alpha_tau_testf3)+0.05, legend=leg, lwd=2, col=colors)
legend(1.6, max(alpha_tau_testf3)-0.5, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)



leg = c(expression(gamma[4][1]), expression(gamma[4][2]), expression(gamma[4][3]), expression(gamma[4][4]), expression(gamma[4][5]), 
        expression(gamma[4][6]), expression(gamma[4][7]), expression(gamma[4][8]), expression(gamma[4][9]))
par(mfrow=c(1,1))
plot(tauv,alpha_tau_testf4[,1],type='l',lwd=2, col=colors[1], xlab=expression(tau), ylab='L2 norm', 
     ylim=c(min(alpha_tau_testf4)-0.05, max(alpha_tau_testf4)+0.05), main=expression(alpha[4](T)))
for(i in 2:(N+1)){
  lines(tauv,alpha_tau_testf4[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(1.8, max(alpha_tau_testf4)+0.05, legend=leg, lwd=2, col=colors)
legend(1.6, max(alpha_tau_testf4)-0.5, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)




leg = c(expression(gamma[5][1]), expression(gamma[5][2]), expression(gamma[5][3]), expression(gamma[5][4]), expression(gamma[5][5]), 
        expression(gamma[5][6]), expression(gamma[5][7]), expression(gamma[5][8]), expression(gamma[5][9]))
par(mfrow=c(1,1))
plot(tauv,alpha_tau_testf5[,1],type='l',lwd=2, col=colors[1], xlab=expression(tau), ylab='L2 norm', 
     ylim=c(min(alpha_tau_testf5)-0.05, max(alpha_tau_testf5)+0.05), main=expression(alpha[5](T)))
for(i in 2:(N+1)){
  lines(tauv,alpha_tau_testf5[,i], type='l',lwd=2, col=colors[i+1])
}
abline(v=tau_TLP, col = "red", lty = 2, lwd=2)
legend(1.8, max(alpha_tau_testf5)+0.05, legend=leg, lwd=2, col=colors)
legend(1.6, max(alpha_tau_testf5)-0.5, legend=c("Selected by BIC"), col=c("red"), lty=2, lwd=2, cex=0.8)

