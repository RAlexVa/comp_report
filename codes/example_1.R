### Source functions
#rm(list=ls())
setwd('..')
source(file.path(getwd(),'functions','IITupdate.R')) #Functions for IIT update
source(file.path(getwd(),'functions','MHupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','MH-IITupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','RN-IITupdate.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','MTM.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','balancing_functions.R')) #Balancing functions
example1 <- function(num_sim=50,
                     max_iter=500*1000,
                     p=500,
                     p_1,
                     theta,
                     threshold=0.1,
                     h,
                     update_step,
                     name_alg){
  simulation_name <- paste0('ex1_',name_alg,'_t',theta,'_p',p,'_p1_',p_1,'_sim',num_sim,'.csv')
  ### Call and define functions to use for the simulation
  X_mode <- c(rep(1,p_1),rep(0,p-p_1)) #Global mode for example 1
  # pi <- function(X){ #Function to return the unnormalized target pi
  #   dif <- -sum(abs(X-X_mode))
  #   return(exp(theta*dif))
  # }
  pi.distribution <- function(X){ #Function to return the log probabilities
    dif <- -sum(abs(X-X_mode))
    return(theta*dif)
  }
  
  # To measure convergence
  pi_F_true <- (choose(p,0:p)*exp(-theta*(0:p)))/((1+exp(-theta))^p) #Vector with the true PI (pushforward measure)
  
  dist_pi <- function(est){ #Fucntion to return the TVD between the true pi and the estimated pi (pushforward measures)
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
      hamm_dist <- sum(abs(X_mode-X)) #distance to mode
      pi_F_est[hamm_dist + 1] <- pi_F_est[hamm_dist + 1] + W #Assign estimated weight
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
  write.csv(cbind(iter_conv,calls_for_pi,last_F_value),file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}
ex1_RN_IIT <- function(num_sim=50,
                     max_iter=500*1000,
                     p=500,
                     p_1,
                     theta,
                     threshold=0.1,
                     h,
                     update_step=RN_IITupdate_log,
                     m=100,
                     name_alg){
  simulation_name <- paste0('ex1_',name_alg,'_t',theta,'_p',p,'_p1_',p_1,'_sim',num_sim,'.csv')
  ### Call and define functions to use for the simulation
  X_mode <- c(rep(1,p_1),rep(0,p-p_1)) #Global mode for example 1

  pi.distribution <- function(X){ #Function to return the log probabilities
    dif <- -sum(abs(X-X_mode))
    return(theta*dif)
  }
  
  # To measure convergence
  pi_F_true <- (choose(p,0:p)*exp(-theta*(0:p)))/((1+exp(-theta))^p) #Vector with the true PI (pushforward measure)
  
  dist_pi <- function(est){ #Fucntion to return the TVD between the true pi and the estimated pi (pushforward measures)
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
    #Initial neighborhood
    neighbors <- sample(1:p,m, replace = F) #Choose m neighbors
    ### Within that for loop need a loop for the steps (considering max number of iterations)
    for(step in 1:max_iter){
      if(step %% 1000 ==0){print(paste(name_alg,'theta',theta,',Iteration:',step,',Simulation:',i))}
      iter <- update_step(X,pi.distribution,h,p,m,neighbors)
      W <- max(iter[[2]],.Machine$double.eps) #Estimated weight of previous state, bounding it from below
      count_PIs <- count_PIs+iter[[3]]
      hamm_dist <- sum(abs(X_mode-X)) #distance to mode
      pi_F_est[hamm_dist + 1] <- pi_F_est[hamm_dist + 1] + W #Assign estimated weight
      
      if(dist_pi(pi_F_est/sum(pi_F_est)) < threshold){ #Compare true dist with normalized est dist.
        #If the chain converged
        iter_conv <- c(iter_conv,step) #Register number of steps needed
        early_finish <- TRUE
        break; #Finish this simulation
      }
      #If the chain has not converged
      X <- iter[[1]] #Update state
      neighbors <- iter[[4]] #Update neighborhood
    }
    calls_for_pi <- c(calls_for_pi,count_PIs) #Register number of times we used the PI function
    last_F_value <- c(last_F_value,dist_pi(pi_F_est/sum(pi_F_est))) #Register last distance from true pi(F)
    #If after max_iter iterations it didn't converge
    if(!early_finish){iter_conv <- c(iter_conv,max_iter)} #register max_iter  
  }
  #iter_conv
  write.csv(cbind(iter_conv,calls_for_pi,last_F_value),file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}

# Parameters
thetas <- c(6,7,8,9) 
p1s <- c(20,50)
set.seed(1348) #Define seed
for(i in 1:length(thetas)){
  for(j in 1:length(p1s)){
    theta_selected <- thetas[i]
    p1_selected <- p1s[j]
    
    example1(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             update_step=IITupdate_log,
             name_alg='IIT')
    
    example1(theta=theta_selected,p_1=p1_selected,
             h=hmin_log,
             update_step=MHupdate_log,
             name_alg='MH')
    
    MH_IIT_defined <- function(X,pi,h,p){
      return(MH_IITupdate_log(X,pi,h,p,rho=0.025))
    }
    example1(theta=theta_selected,p_1=p1_selected,
             h=hmin_log,
             update_step=MH_IIT_defined,
             name_alg='MH-IIT')
    
    ex1_RN_IIT(theta=theta_selected,p_1=p1_selected,
             h=hmin_log,
             max_iter = 500*1000,
             name_alg='RN-IIT')
    
    MTM_defined <- function(X,pi,h,p){
      return(MTMupdate_log(X,pi,h,p,m=100))
    }
    example1(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             update_step=MTM_defined,
             name_alg='MTM')
  }
}




