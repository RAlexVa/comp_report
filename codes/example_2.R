#rm(list=ls())
setwd('..')
### Source functions
source(file.path(getwd(),'functions','IITupdate.R')) #Functions for IIT update
source(file.path(getwd(),'functions','MHupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','MH-IITupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','RN-IITupdate.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','MTM.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','balancing_functions.R')) #Balancing functions
example2 <- function(num_sim=50,
                     max_iter=300*1000,
                     p=500,
                     theta,
                     threshold=0.2,
                     h,
                     update_step,
                     name_alg){
  simulation_name <- paste0('ex2_',name_alg,'_t',theta,'_p',p,'_sim',num_sim,'.csv')
  ### Call and define functions to use for the simulation
  X_mode <- c(1,rep(0,p-1)) #Global mode for example 2
  pi.distribution <- function(X){ #Function to return the log probabilities
    if(X[1]==1){ell <- (sum(X)-1)}
    if(X[1]==0){ell <- 2*p-sum(X)}
    return(-theta*ell)
  }
  Fm <- function(X){
    if(X[1]==1){val <- (sum(X)-1)}
    if(X[1]==0){val <- p}
    return(val)
  }
  # To measure convergence
  pi_F_true <- (choose(p-1,0:(p-1))*exp(-theta*(0:(p-1)))) #Vector with the true PI (pushforward measure)
  pi_F_true <- c(pi_F_true,exp(log(1+exp(-theta))*(p-1)-theta*(p+1))) #Adding the last coordinate
  pi_F_true <- pi_F_true/(exp(log(1+exp(-theta))*(p-1)-theta*(p+1)) + exp((p-1)*log(1+exp(-theta))))
  
  
  dist_pi <- function(est){ #Function to return the TVD between the true pi and the estimated pi (pushforward measures)
    if(length(est)!=length(pi_F_true)){
      print('Size of vector is not the same')
      return(1)}
    return(sum(abs(est - pi_F_true)))
  }
  
  ### Need a for loop considering the number of simulations
  iter_conv <- c()
  calls_for_pi <- c()
  last_F_value <- c()
  for(i in 1:num_sim){
    print(paste(name_alg,'Simulation:',i,'theta',theta))
    ### To check thresholds
    pi_F_est <- numeric(p+1)
    early_finish <- FALSE #To check if convergence was achieved before max. iterations
    #Initialize
    X <- rep(0,p)
    #Counting the number of times we used the PI function
    count_PIs <- 0
    ### Within that for loop need a loop for the steps (considering max number of iterations)
    for(step in 1:max_iter){
      if(step %% 1000 ==0){print(paste(name_alg,'theta',theta,',Iteration:',step,',Simulation:',i))}
      iter <- update_step(X,pi.distribution,h,p)
      W <- max(iter[[2]],.Machine$double.eps) #Estimated weight of previous state, bounding it from below
      count_PIs <- count_PIs+iter[[3]]
      #Adding a minimum that is above 0 to avoid issues with very low numbers
      Feval <- Fm(X) +1 #Function F evaluated
      pi_F_est[Feval] <- pi_F_est[Feval] + W #Assign estimated weight

      if(dist_pi(pi_F_est/sum(pi_F_est)) < threshold){ #Compare true dist with normalized est dist.
        #If the chain converged
        iter_conv <- c(iter_conv,step) #Register number of steps needed
        early_finish <- TRUE
        break; #Finish this simulation
      }
      #If the chain has not converged
      X <- iter[[1]] #Update state
    }
    calls_for_pi <- c(calls_for_pi,count_PIs) #Register number of times we used the PI function
    last_F_value <- c(last_F_value,dist_pi(pi_F_est/sum(pi_F_est))) #Register last distance from true pi(F)
    #If after max_iter iterations it didn't converge
    if(!early_finish){iter_conv <- c(iter_conv,max_iter)} #register max_iter  
  }
  #iter_conv
  write.csv(cbind(iter_conv,calls_for_pi,last_F_value),file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}

MTM_defined <- function(X,pi,h,p){
  return(MTMupdate_log(X,pi,h,p,m=100))
}

# Parameters
thetas <- c(6,7,8,9) 
set.seed(348) #Define seed
for(i in 1:length(thetas)){
    theta_selected <- thetas[i]
    
    example2(theta=theta_selected,
             h=hsq_log,
             update_step=IITupdate_log,
             name_alg='IIT')
    
    example2(theta=theta_selected,
             h=hmin_log,
             update_step=MHupdate_log,
             name_alg='MH')
    
    MH_IIT_defined <- function(X,pi,h,p){
      return(MH_IITupdate_log(X,pi,h,p,rho=0.025))
    }
    example2(theta=theta_selected,
             h=hmin_log,
             update_step=MH_IIT_defined,
             name_alg='MH-IIT')
    
    RN_IIT_defined <- function(X,pi,h,p){
      return(RN_IITupdate_log(X,pi,h,p,m=100))
    }
    example2(theta=theta_selected,
             h=hmin_log,
             update_step=RN_IIT_defined,
             name_alg='RN-IIT')
    
    MTM_defined <- function(X,pi,h,p){
      return(MTMupdate_log(X,pi,h,p,m=100))
    }
    example2(theta=theta_selected,
             h=hsq_log,
             update_step=MTM_defined,
             name_alg='MTM')
  }
