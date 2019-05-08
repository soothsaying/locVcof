

####X is the covariate
softth<-function(b,lam)
{
absb=abs(b)
if(absb>lam) b1=absb-lam
else b1=0
if(b>=0) return(b1)
else return(-b1)
}

####Y is the response variable
####N is the size of group variables
groupLASSO<-function(x,y,N,lambda,weight)
{
d<-dim(x)[2]
n<-length(y)
p<-d/N

#### First standadize the data
meanx=apply(x,2,mean)
meany=mean(y)
y=y-meany
for(i in 1:d)
x[,i]=x[,i]-meanx[i]

Transform<-matrix(0,N,d)
for(i in 1:p)
{
temp=((i-1)*N+1):(i*N)
sigma<-(t(x[,temp])%*%x[,temp])
e<-eigen(sigma)
v<-e$vectors
sigma_half<-v %*% diag(sqrt(e$values)) %*% t(v) 
x[,temp]=x[,temp]%*%solve(sigma_half)
Transform[,temp]=solve(sigma_half)
}

beta_ols<-solve(t(x)%*%x)%*%t(x)%*%y

diff=1
epsilon=10^(-6)
beta_Lasso=beta_ols

while(diff>=epsilon)
{
temp_beta=beta_Lasso
for(i in 1:p)
{
temp=((i-1)*N+1):(i*N)
y_i=y-x[,-temp]%*%beta_Lasso[-temp]
x_i=x[,temp]
b_i=solve(t(x_i)%*%x_i)%*%t(x_i)%*%y_i
normb_i=mean((x_i%*%b_i)^2)
beta_Lasso[temp]=softth(normb_i,lambda*weight[i])*b_i/normb_i
}
diff=max(abs(beta_Lasso-temp_beta))
}

for(i in 1:p)
{
temp=((i-1)*N+1):(i*N)
beta_Lasso[temp]=Transform[,temp]%*%beta_Lasso[temp]
}

return(beta_Lasso)
}