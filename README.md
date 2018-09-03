## Bayesian-Single-Index-Model

# BSIM without missing covariates
 Codes_for_bsim_no_miss.R contains codes when none of the covariates are missing. 
 These codes are written when the number of covariates is p=3. 
 Example is given in example_bsim_no_miss.R file. 

 bsim_no_miss function has 5 inputs and they are described as follows. 
# Inputs
 x Covariates p=3
 
 y Response
 
 N Number of MCMC iterations. 
 
 burn Number of burn-ins+1
 
 grid.width grid width used for the grid search. A smaller value will make computation more accurate but more time consuming. 
 
# Outputs
 theta_est Estimated credible intervals of two theta values 
 storetheta Values of theta over N mcmc iterations. 


# BSIM with missing covariates

 Codes_for_bsim_with_miss.R contains codes when there are missing observations in the first covariate.  
 These codes are written when the number of covariates is p=3. 
 Example is given in example_bsim_with_missing.R file. 

bsim_miss_cov function has 5 inputs and they are described as follows. 
#Inputs
x Covariates p=3
y Response
N Number of MCMC iterations. 
burn  Number of burn-ins+1
grid.width  grid width used for the grid search. A smaller value will make computation more accurate but more time consuming. 
#Outputs
theta_est Estimated 0.025, 0.5 and 0.975 percentiles of the two thetas based on complete part of the data.
thetaci Estimated 0.025, 0.5 and 0.975 percentiles after incorporating missing covariates.
storetheta  Values of theta over N mcmc iterations. 
weight  Contains weights of corresponding to storetheta values. 



