rm(list=ls())

library(MASS)
library(Matrix)
library(mvtnorm)
library(foreach)
library(doSNOW)
library(coda)

n_all=100    # Sample Size
p=3		#dimension of X

mu=c(0,0,0)
rho=0
Cx<-rbind(c(1,rho,rho), c(rho,1,rho), c(rho, rho,1))
x<-mvrnorm(n=n_all,mu=mu,Sigma=Cx,tol=1e-8) #simulation of f without noise
x<-t(x)

alpha<-c(1,1,1)  ####Final Value of alpha 
alpha<-alpha/sqrt(sum(alpha^2))
z<-matrix(0,n_all) #z=alpha'x
z<-t(x)%*%alpha



sigma<-.3
f<-exp(z)
y<-f + rnorm(n_all,0,sd=sigma) # adding noise
# Data Generated, (y,x) the dataset at this point

y<-y-mean(y)


var(f)

f_all<-f
x_all<-x
z_all<-z
y_all<-y
data_frame<-cbind(t(x_all), y, f)
colnames(data_frame)=c('x1', 'x2', 'x3', 'y','f')                      
xobs<-x
yobs<-y
n=n_all
p=3		#dimension of X



result<-bsim_no_miss(x=(xobs), y=yobs, N=50, burn=21, grid.width=.1)
