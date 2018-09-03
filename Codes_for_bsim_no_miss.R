###Code for BSIM


rm(list=ls())


library(MASS)
library(Matrix)
library(mvtnorm)
library(foreach)
library(doSNOW)
library(coda)

assign("last.warning", NULL, envir = baseenv())
#set.seed(5000)

##Initializing for parallel computing
#cores<-20
#cl <- makeCluster(cores, type = "SOCK")

#registerDoSNOW(cl)


###The function that needs to be maximized at grid points
optimized<-function(f,xobs,kappa,thetahat,p)       
{
  n=length(f)
  thetahat<-as.numeric(as.matrix(thetahat))
  alphahat<-matrix(0,p,1)
  alphahat[1]<-sin(thetahat[1])
  alphahat[2]<-sin(thetahat[2])*cos(thetahat[1])
  alphahat[3]<-cos(thetahat[2])*cos(thetahat[1])
  #alphahat<-as.matrix(alphahat)
  x<-t(xobs)%*%(alphahat)
  f<-f[order(x)]
  x<-x[order(x)]
  
  r=matrix(0,n+1,1)
  r[2:n]=exp(-kappa*(x[2:n]-x[1:(n-1)]))
  r[1]=0
  r[n+1]=0
  
  e<-matrix(0,n+1,1)
  e<-r/(1-(r^2))
  d<-matrix(0,n,1)
  d[1:n]<-1+r[1:n]*e[1:n]+r[2:(n+1)]*e[2:(n+1)]
  
  Cinv<-matrix(0,n,n)
  diag(Cinv)=d
  
  for(i in 1:(n-1))
  {
    Cinv[i,i+1]=-e[i+1]
    Cinv[i+1,i]=-e[i+1]
  }
  
  
  return(0.5*log((det(Cinv)))+(-0.5*f%*%Cinv%*%f))
  
  
}

#The function used to find the parameters of Beta distribution. 
#Given the values of \theta

betaalpha1<-function(theta)
{
  p<-length(theta)
  t<-matrix(0,p,1)
  t[1]<-theta[1]/(pi)
  t[2]<-theta[2]/pi
  abeta<-matrix(0,p,1) #Higher alpha is lower variance
  abeta[1]<-5000 
  abeta[2]<-5000 
  bbeta<-matrix(0,p,1)
  for(i in 1:p)
  {bbeta[i]<-(abeta[i]-1-abeta[i]*t[i]+2*t[i])/t[i]}
  return(cbind(abeta,bbeta))
}
thetatoalpha<-function(theta)
{
  alpha<-matrix(0, length(theta)+1)
  alpha[1]<-sin(theta[1])
  alpha[2]<-(cos(theta[1]))*sin(theta[2])
  alpha[3]<-(cos(theta[1]))*cos(theta[2])
  return(alpha)
}

alphatotheta<-function(alpha)
{
  theta<-matrix(0, length(alpha)-1)
  theta[1]<-asin(alpha[1])
  theta[2]<-asin(alpha[2]/(cos(theta[1])))
}



###Final Code

bsim_no_miss<-function(x, y, N, burn, grid.width=.5)
{
  tm<-proc.time() #Finding duration 
  
  xobs<-x
  yobs<-y
  
  n<-dim(xobs)[2]
  p<-dim(xobs)[1]
  
  theta<-matrix(0,p-1,1) #Matrix containing theta's
  
  theta[1]<-runif(1,0,pi) #Initial distribution
  theta[2]<-runif(1,0,pi)
  thetat<-theta
  parabeta<-betaalpha1(thetat)
  
  
  #Converting theta to alpha
  alphahat<-matrix(0,p,1)
  alphahat[1]<-sin(theta[1])
  alphahat[2]<-sin(theta[2])*cos(theta[1])
  alphahat[3]<-cos(theta[2])*cos(theta[1])
  
  x<-t(xobs)%*%alphahat #x'\alpha
  storeorder<-order(x) 
  #Ordering x and y for fitting OU process
  y<-yobs[order(x)]  
  x<-x[order(x)]
  
  kappa=2 # scale parameter for C matrix, Initialization
  
  #Finding inverse of C
  r=matrix(0,n+1,1)
  r[2:n]=exp(-kappa*(x[2:n]-x[1:(n-1)]))
  r[1]=0
  r[n+1]=0
  
  e<-matrix(0,n+1,1)
  e<-r/(1-(r^2))
  d<-matrix(0,n,1)
  d[1:n]<-1+r[1:n]*e[1:n]+r[2:(n+1)]*e[2:(n+1)]
  
  Cinv<-matrix(0,n,n)
  diag(Cinv)=d
  
  for(i in 1:(n-1))
  {
    Cinv[i,i+1]=-e[i+1]
    Cinv[i+1,i]=-e[i+1]
  }
  
  
  #generate from f
  sigma2<-(0.01^2) #Initial sigma^2
  Sigma<-solve(Cinv+diag(n)/sigma2)#+diag(n)*0.005
  mu<-(Sigma%*%y)/sigma2
  
  
  f<-mvrnorm(1,mu=mu,Sigma=Sigma,tol=1e-8) # Better to Use method from Rue
  
  #Generate from sigma2
  a=1;
  b=0.01;
  err<-sum((y-f)^2)/2
  sigma2<-1/rgamma(1,shape=(a+n/2),rate=(b+err))
  sigma2<-sigma^2
  #selecting kappa for getting proper C matrix
  kappa=2
  
  
  kappas<-as.matrix(seq(0.5,4,0.05))
  probs<-matrix(0,length(kappas))
  
  T<-(x%*%t(matrix(1,n,1)))
  
  a=1;
  b=0.01;
  thetas<-as.matrix(expand.grid(x=seq(0.01,pi-0.01,grid.width),y=seq(0.01,pi-0.01,grid.width)))
  num.calc<-dim(thetas)[1]
  g<-f[order(storeorder)] 
  
  ##Range of theta is defined
  PP1<-apply(thetas,1,function(x) optimized(g,xobs,kappa,x,p)) #Evaluation of function for elements in PP1
  theta=thetas[which.max(PP1),] #Finding the maximum 
  
  thetat<-theta
  parabeta<-betaalpha1(thetat)
  
  
  #Converting theta to alpha
  alphahat<-matrix(0,p,1)
  alphahat[1]<-sin(theta[1])
  alphahat[2]<-sin(theta[2])*cos(theta[1])
  alphahat[3]<-cos(theta[2])*cos(theta[1])
  
  x<-t(xobs)%*%alphahat #x'\alpha
  storeorder<-order(x) 
  #Ordering x and y for fitting OU process
  y<-yobs[order(x)]  #### Check  plot(xobs[1,],y[order(storeorder)])= plot(xobs[1,],yobs)
  x<-x[order(x)]
  
  kappa=2 # scale parameter for C matrix, Initialization
  
  #Finding inverse of C
  r=matrix(0,n+1,1)
  r[2:n]=exp(-kappa*(x[2:n]-x[1:(n-1)]))
  r[1]=0
  r[n+1]=0
  
  e<-matrix(0,n+1,1)
  e<-r/(1-(r^2))
  d<-matrix(0,n,1)
  d[1:n]<-1+r[1:n]*e[1:n]+r[2:(n+1)]*e[2:(n+1)]
  
  Cinv<-matrix(0,n,n)
  diag(Cinv)=d
  
  for(i in 1:(n-1))
  {
    Cinv[i,i+1]=-e[i+1]
    Cinv[i+1,i]=-e[i+1]
  }
  
  
  #generate from f
  sigma2<-(0.01^2) #Initial sigma^2
  Sigma<-solve(Cinv+diag(n)/sigma2)#+diag(n)*0.005
  mu<-(Sigma%*%y)/sigma2
  
  f<-mvrnorm(1,mu=mu,Sigma=Sigma,tol=1e-8) # Better to Use method from Rue
  
  #Generate from sigma2
  a=7;
  b=2;
  err<-sum((y-f)^2)/2
  sigma2<-1/rgamma(1,shape=(a+n/2),rate=(b+err))
  sigma2<-sigma^2
  
  kappa=2
  
  
  kappas<-as.matrix(seq(0.5,4,0.05))
  probs<-matrix(0,length(kappas))
  T<-(x%*%t(matrix(1,n,1)))
  
  a=2;
  b=0.01;
  
  N=N  ##Total number of MCMC iterations
  
  counter=0
  parabeta1<-matrix(0,2,1)
  storef<-matrix(0,n,N)
  storetheta<-matrix(0,2,N)
  storekappa<-matrix(0,N)
  storeprop<-matrix(0,N)
  storecounter<-matrix(0,N)
  storesigma2<-matrix(0,N)
  avg<-matrix(0,n)
  for ( rep in 1:N)
  {
    
    
    r=matrix(0,n+1,1)
    r[2:n]=exp(-kappa*(x[2:n]-x[1:(n-1)]))
    r[1]=0
    r[n+1]=0
    
    e<-matrix(0,n+1,1)
    e<-r/(1-(r^2))
    d<-matrix(0,n,1)
    d[1:n]<-1+r[1:n]*e[1:n]+r[2:(n+1)]*e[2:(n+1)]
    
    Cinv<-matrix(0,n,n)
    diag(Cinv)=d
    
    for(i in 1:(n-1))
    {
      Cinv[i,i+1]=-e[i+1]
      Cinv[i+1,i]=-e[i+1]
    }
    
    #generate from f
    
    Sigma<-solve(Cinv+diag(n)/sigma2)#+diag(n)*0.005
    mu<-(Sigma%*%y)/sigma2
    f<-mvrnorm(1,mu=mu,Sigma=Sigma,tol=1e-8) # Better to Use method from Rue
    
    #Generate from sigma2
    err<-sum((y-f)^2)/2
    sigma2<-1/rgamma(1,shape=(a+n/2),rate=(b+err))
    storesigma2[rep]<-sigma2
    
    #selecting kappa for getting proper C matrix
    ep<-matrix(0,length(kappas))
    for (pt in 1:length(kappas))
    {
      kappa=kappas[pt]
      #C<-kappa*abs(T-t(T))
      #C<-exp(-C)
      #probs[p]<-dmvnorm(f,mean=matrix(0,n,1),sigma=C,log=TRUE)
      r=matrix(0,n+1,1)
      r[2:n]=exp(-kappa*(x[2:n]-x[1:(n-1)]))
      r[1]=0
      r[n+1]=0
      
      e<-matrix(0,n+1,1)
      e<-r/(1-(r^2))
      d<-matrix(0,n,1)
      d[1:n]<-1+r[1:n]*e[1:n]+r[2:(n+1)]*e[2:(n+1)]
      
      Cinv<-matrix(0,n,n)
      diag(Cinv)=d
      
      for(i in 1:(n-1))
      {
        Cinv[i,i+1]=-e[i+1]
        Cinv[i+1,i]=-e[i+1]
      }
      probs[pt]<-(0.5*log(det(Cinv))-0.5*f%*%Cinv%*%f)
    }
    
    probs<-probs-max(probs)
    probs<-exp(probs)
    
    kas<-which(rmultinom(1,1,probs)==1)   #Selecting kappa from Multinomial setup
    kappa<-kappas[kas]
    storekappa[rep]<-kappa
    
    
    g<-f[order(storeorder)] 
    
    ##Range of theta is defined
    #PP1<-apply(thetas,1,function(x) optimized(g,xobs,kappa,x,p)) #Evaluation of function for elements in PP1
    #theta=thetas[which.max(PP1),] #Finding the maximum 
    idxx<-((rep%%10)/100+0.01)
    thetas<-as.matrix(expand.grid(x=seq(idxx,pi-0.01,grid.width),y=seq(idxx,pi-0.01,grid.width)))
    num.calc<-dim(thetas)[1]
    PP1<-apply(thetas,1,function(x) optimized(g,xobs,kappa,x,p)) #Evaluation of function for elements in PP1
    theta=thetas[which.max(PP1),] #Finding the maximum 
    
    #calc<-foreach(i=1:num.calc, .combine="rbind") %dopar%{
    #optimized(g,xobs,kappa,thetas[i,],p)}
    
    #theta=thetas[which.max(calc),]
    #thetat is the old theta,  theta is mode of new theta
    
    ####MH algorithm#####
    Stheta<-0.1*diag(1,n)
    g<-mvrnorm(1,mu=f,Sigma=Stheta)
    ##Need to find Rtheta=Stheta-Stheta%*%solve(Stheta+solve(invC))%*%Stheta
    ##xobs and theta are given. 
    
    alpha[1]<-sin(thetat[1])
    alpha[2]<-sin(thetat[2])*cos(thetat[1])
    alpha[3]<-cos(thetat[2])*cos(thetat[2])
    alpha<-alpha/sqrt(sum(alpha^2))
    z<-matrix(0,n) #z=alpha'x
    for (i in 1:n)
    {
      z[i]<-xobs[,i]%*%alpha      #these z's are not ordered as g, i.e., as f
    }
    
    z<-z[(storeorder)]
    R1<-(matrix(1,n))%*%t(z)
    R<-(abs(R1-t(R1)))
    C<-exp(-kappa*R)
    
    Sigmatheta<-C
    
    Rtheta<-Stheta-Stheta%*%solve(Stheta+Sigmatheta)%*%Stheta
    
    LRtheta<-t(chol(Rtheta))
    mtheta_g<-Rtheta%*%solve(Stheta)%*%g
    
    eta<-solve(LRtheta)%*%(f-mtheta_g)
    
    ##Generating proposal theta, named as theta1
    parabeta1<-betaalpha1(theta)        #theta is the point obtained from maximizing over grids
    theta1<-matrix(0,2,1) #New possible theta
    theta1[1]<-pi*rbeta(1,parabeta1[1,1],parabeta1[1,2])
    theta1[2]<-pi*rbeta(1,parabeta1[2,1],parabeta1[2,2])
    
    ##Finding the values of LRthetaprime, mtheta_gprime for new value of theta, theta1
    alpha1<-matrix(0,3)
    alpha1[1]<-sin(theta1[1])
    alpha1[2]<-sin(theta1[2])*cos(theta1[1])
    alpha1[3]<-cos(theta1[2])*cos(theta1[2])
    z1<-matrix(0,n) #z=alpha'x
    for (i in 1:n)
    {
      z1[i]<-xobs[,i]%*%alpha1
    }
    z1<-z1[(storeorder)]
    
    R1<-(matrix(1,n))%*%t(z1)
    R<-(abs(R1-t(R1)))
    C<-exp(-kappa*R)
    
    Sigmathetaprime<-C
    
    Rthetaprime<-Stheta-Stheta%*%solve(Stheta+Sigmathetaprime)%*%Stheta
    
    LRthetaprime<-t(chol(Rthetaprime))
    mtheta_gprime<-Rtheta%*%solve(Stheta)%*%g
    
    fprime<-LRthetaprime%*%eta+mtheta_gprime
    
    
    
    
    L<-(dmvnorm(as.numeric(fprime),mean=rep(0,n),sigma=Sigmathetaprime,log=TRUE)
        - dmvnorm(as.numeric(f), mean=rep(0,n), sigma=Sigmatheta, log=TRUE)
        + dmvnorm(as.numeric(g), mean=rep(0,n), sigma=(Sigmathetaprime+Rthetaprime), log=TRUE)
        - dmvnorm(as.numeric(g),mean=rep(0,n), sigma=(Sigmatheta+Rtheta),log=TRUE)
        - dbeta((theta1[1]/pi), parabeta[1,1], parabeta[1,2], log=TRUE)
        - dbeta((theta1[2]/pi), parabeta[2,1], parabeta[2,2], log=TRUE)
        + dbeta((thetat[1]/pi), parabeta1[1,1], parabeta1[1,2], log=TRUE)
        + dbeta((thetat[2]/pi), parabeta1[2,1], parabeta1[2,2], log=TRUE))
    + (2*log(abs(cos(theta1[1])))) - log(abs(cos(theta1[2])))
    - (2*log(abs(cos(thetat[1])))) - log(abs(cos(thetat[2])))
    ratio<-exp(L)
    sim<-runif(1)
    if(sim<ratio)
    {
      theta<-theta1; 
      thetat<-theta1;
      f<-fprime
      storecounter[rep]=storecounter[rep]+1;
      parabeta<-parabeta1;
    }
    if(sim >ratio){theta<-thetat}
    
    storetheta[,rep]<-theta
    storef[storeorder,rep]<-f
    print(rep)
    
    
    
    alphahat<-matrix(0,p,1)
    alphahat[1]<-sin(theta[1])
    alphahat[2]<-sin(theta[2])*cos(theta[1])
    alphahat[3]<-cos(theta[2])*cos(theta[1])
    
    
    x<-t(xobs)%*%alphahat  	#Updating using in alphahat
    storeorder<-order(x)	#Storing the ordering
    y<-yobs[order(x)]		#Updating ordered y
    x<-x[order(x)]		#Updating ordered x
    #print(kappa)
    #print(theta)
    print(proc.time()-tm)
  }
  
  avg<-avg
  burn=burn   ###Burn-in
  # s<-summary(mcmc(t(storetheta[,burn:N])))
  # capture.output(s,file="exp_100_0_5.txt",append=TRUE)
  # s<-rejectionRate(mcmc(t(storetheta[,burn:N])))
  #  capture.output(s,file="exp_100_0_5.txt",append=TRUE)
  
  # write.table("\n \n \n ", file="exp_100_0_5.txt",append=TRUE)
  
  
  
  #stopCluster(cl) #Closing Cluster
  completedatathetaci<-matrix(0, 2,3)
  completedatathetaci[1,]=as.matrix(quantile(storetheta[1,burn:N], c(0.025, .5,.975)))
  completedatathetaci[2,]=as.matrix(quantile(storetheta[2,burn:N], c(0.025, .5,.975)))
  theta_est<-rbind(completedatathetaci[1,], completedatathetaci[2,])
  
  coverprob1<- ifelse((quantile(storetheta[1,burn:N],0.025)<0.6154797) & (quantile(storetheta[1,burn:N],.975)>0.6154797),1,0)
  coverprob2<- ifelse((quantile(storetheta[2,burn:N],0.025)<0.7853982) & (quantile(storetheta[2,burn:N],.975)>0.7853982),1,0)
  
  meantheta1<-mean(storetheta[1,burn:N])
  meantheta2<-mean(storetheta[2,burn:N])
  
  
  mediantheta1<-median(storetheta[1,burn:N])
  mediantheta2<-median(storetheta[2,burn:N])
  
  return(theta=theta_est, storetheta=storetheta)
  
  
  
}
