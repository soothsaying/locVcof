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

SCADth<-function(b,lambda,a)
{
absb=abs(b)
c1=((a-1)*absb-a*lambda)/((a-2)*absb)
c2=(absb-lambda)*(absb>=lambda)/absb
if(absb>=a*lambda)
b1=b
else if ((absb<a*lambda)&(absb>=2*lambda))
b1=c1*b
else
b1=c2*b
return(b1)
}


####Y is the response variable
####N is the size of group variables
groupSCAD<-function(x,y,N,lambda,a)
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
beta_SCAD=beta_ols

step=0
while((diff>=epsilon)&(step<=100))
{
temp_beta=beta_SCAD
for(i in 1:p)
{
temp=((i-1)*N+1):(i*N)
y_i=y-x[,-temp]%*%beta_SCAD[-temp]
x_i=x[,temp]
b_i=solve(t(x_i)%*%x_i)%*%t(x_i)%*%y_i
normb_i=mean((x_i%*%b_i)^2)
normb_i2=SCADth(normb_i,lambda,a)/normb_i
beta_SCAD[temp]=normb_i2*b_i
}
diff=mean((beta_SCAD-temp_beta)^2)
step=step+1
}

if(step>100) print(step)

for(i in 1:p)
{
temp=((i-1)*N+1):(i*N)
beta_SCAD[temp]=Transform[,temp]%*%beta_SCAD[temp]
}
return(beta_SCAD)

}
