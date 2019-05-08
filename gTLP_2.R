#library("lars")
####X is the covariate
softth<-function(b,lam)
{
absb=abs(b)
if(absb>lam) b1=absb-lam
else b1=0
if(b>=0) return(b1)
else return(-b1)
}

TLPth<-function(b,lambda,tau)
{
absb=abs(b)
r=lambda/tau
if((absb>=max(r/2+tau,r))||((r>=absb)&(absb>=max(sqrt(2*lambda),tau))))
b1=b
else if ((absb<tau)||((absb>=r)&(absb<=(r/2+tau))))
b1=softth(b,r)
else
b1=0
return(b1)
}

####Y is the response variable
####N is the size of group variables
groupTLP<-function(x,y,N,lambda,tau)
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

#lasso.res<-lars(x,y, type="lasso", use.Gram=F)
#beta_ols<-solve(t(x)%*%x)%*%t(x)%*%y
beta_ols<-rep(0,d)


diff=1
epsilon=10^(-6)
beta_TLP=beta_ols

step=0
while(diff>=epsilon & step<=100)
{
temp_beta=beta_TLP
ind=0
for(i in 1:p)
{
temp=((i-1)*N+1):(i*N)
y_i=y-x[,-temp]%*%beta_TLP[-temp]
x_i=x[,temp]
b_i=solve(t(x_i)%*%x_i)%*%t(x_i)%*%y_i
normb_i=sqrt(mean((x_i%*%b_i)^2))
normb_i2=TLPth(normb_i,lambda,tau)/normb_i
beta_TLP[temp]=normb_i2*b_i
#print(normb_i)
#print(normb_i2)
}
diff=mean((beta_TLP-temp_beta)^2)
step=step+1
#print(step)
#print(diff)
}


if(step>100) print(step)

for(i in 1:p)
{
temp=((i-1)*N+1):(i*N)
beta_TLP[temp]=Transform[,temp]%*%beta_TLP[temp]
}
return(beta_TLP)

}
