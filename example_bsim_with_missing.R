library(MASS)
library(Matrix)
library(mvtnorm)
library(foreach)
library(doSNOW)
library(coda)


###Example for Missing covariates

n_all=100    # Sample Size
p=3		#dimension of X
#x<-rnorm(p*n_all,0,1) #x, the covariates
mu=c(0,0,0)
rho=0
Cx<-rbind(c(1,rho,rho), c(rho,1,rho), c(rho, rho,1))
x<-mvrnorm(n=n_all,mu=mu,Sigma=Cx,tol=1e-8) 
x<-t(x)

alpha<-c(1,1,1)  ###alpha_0
alpha<-alpha/sqrt(sum(alpha^2))
z<-matrix(0,n_all) #z=alpha'x
z<-t(x)%*%alpha

sigma<-.5  ##Error Standard deviation

f<-exp(z)  ##True link function
y<-f + rnorm(n_all,0,sd=sigma) # adding noise
# Data Generated, (y,x) the dataset at this point

y<-y-mean(y)
var(f)

f_all<-f
x_all<-x
z_all<-z
y_all<-y

##Missing Data
probs<-1/(1+exp(4+2*x[2,]+3*x[3,]+y))  ###5 as beta_0 gives 10% miss ###4 as beta_0 gives a bit more than 10% miss
delta<-rbinom(n_all, 1, probs) #delta indicates which of the x1's are missing. 
sum(delta)
x[1,(delta==1)]<-NA

xobs<-x

##Fitting data with missing covariates 
result<-bsim_miss_cov(x=t(xobs), y=y_all, N=50, burn=21)
