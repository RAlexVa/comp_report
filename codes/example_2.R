#rm(list=ls())
setwd('..')
### Source functions
source(file.path(getwd(),'functions','IITupdate.R')) #Functions for IIT update
source(file.path(getwd(),'functions','MHupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','MH-IITupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','RN-IITupdate.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','MTM.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','balancing_functions.R')) #Balancing functions
source(file.path(getwd(),'functions','IIT-RFupdate.R')) #Functions for IIT-RF update
example2 <- function(num_sim=50,
                     max_iter=500*1000,
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
    X <- c(rep(0,10),rep(1,p-10))
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
ex2_RN_IIT <- function(num_sim=50,
                     max_iter=500*1000,
                     p=500,
                     theta,
                     threshold=0.2,
                     h,
                     update_step=RN_IITupdate_log,
                     m=100,
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
    X <- c(rep(0,10),rep(1,p-10))
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
ex2_RF_IIT <- function(num_sim=50,
                     max_iter=500*1000,
                     p=500,
                     theta,
                     threshold=0.2,
                     h,
                     update_step=IIT_RFupdate_log,
                     name_alg,
                     initial_K){
  simulation_name <- paste0('ex2_',name_alg,'_t',theta,'_p',p,'_sim',num_sim,'_K',round(initial_K,2),'.csv')
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
  iteration_K <- c() #Specific for RF-IIT
  final_K <- c()#Specific for RF-IIT
  for(i in 1:num_sim){
    print(paste(name_alg,'Simulation:',i,'theta',theta))
    ### To check thresholds
    pi_F_est <- numeric(p+1)
    early_finish <- FALSE #To check if convergence was achieved before max. iterations
    bounding_K <- initial_K #Specific for RF-IIT
    #Initialize
    X <- c(rep(0,10),rep(1,p-10))
    #Counting the number of times we used the PI function
    count_PIs <- 0
    #Counting the number of times we updated K
    count_K <- 0 #Specific for RF-IIT
    ### Within that for loop need a loop for the steps (considering max number of iterations)
    for(step in 1:max_iter){
      if(step %% 1000 ==0){print(paste(name_alg,'theta',theta,',Iteration:',step,',Simulation:',i))}
      iter <- update_step(X,pi.distribution,h,p,bounding_K)#Specific for RF-IIT
      W <- max(iter[[2]],.Machine$double.eps) #Estimated weight of previous state, bounding it from below
      count_PIs <- count_PIs+iter[[3]]
      bounding_K <- iter[[4]]#Specific for RF-IIT
      count_K <- count_K+iter[[5]]#Specific for RF-IIT
      #Adding a minimum that is above 0 to avoid issues with very low numbers
      Feval <- Fm(X) +1 #Function F evaluated
      pi_F_est[Feval] <- pi_F_est[Feval] + W #Assign estimated weight
      if(pi_F_est[Feval]==Inf){pi_F_est[Feval] <- .Machine$double.xmax} #Fix huge values
      fix_dist <- function(pi_F_est){
        #First fix
        totalizer <- sum(pi_F_est)
        if(totalizer==Inf){totalizer <- .Machine$double.xmax}
        newpi <- pi_F_est/totalizer
        #recurrent fix
        while(sum(newpi)==Inf){
          totalizer <- sum(newpi)
          if(totalizer==Inf){totalizer <- .Machine$double.xmax}
          newpi <- newpi/totalizer
        }
        return(newpi)
      }
      pi_F_est <- fix_dist(pi_F_est)

      if(dist_pi(pi_F_est/totalizer) < threshold){ #Compare true dist with normalized est dist.
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
    iteration_K <- c(iteration_K,count_K)#Specific for RF-IIT
    final_K <- c(final_K,bounding_K)#Specific for RF-IIT
    #If after max_iter iterations it didn't converge
    if(!early_finish){iter_conv <- c(iter_conv,max_iter)} #register max_iter  
  }
  #iter_conv
  write.csv(cbind(iter_conv,calls_for_pi,last_F_value),file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}

# MH_IIT_defined <- function(X,pi,h,p){
#   return(MH_IITupdate_log(X,pi,h,p,rho=0.025))
# }
# theta_selected <- 6
# h_defined <- function(r){return(hc_log(r,c=2*theta_selected))}
# 
# num_sim=1;
# max_iter=10*1000;
# p=500;
# theta=6;
# threshold=0.2;
# h=h_defined;
# update_step=MH_IIT_defined;
# name_alg='MH_IIT'



# Parameters
thetas <- c(6,7,8,9) 
set.seed(1235) #Define seed
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
             name_alg='MH-IIT-1')
    
    h_defined <- function(r){return(hc_log(r,c=2*theta_selected))}
    example2(theta=theta_selected,
             h=h_defined,
             update_step=MH_IIT_defined,
             name_alg='MH-IIT-2')
    
    ex2_RN_IIT(theta=theta_selected,
             h=hmin_log,
             update_step=RN_IITupdate_log,
             name_alg='RN-IIT')
    
    MTM_defined <- function(X,pi,h,p){
      return(MTMupdate_log(X,pi,h,p,m=100))
    }
    example2(theta=theta_selected,
             h=hsq_log,
             update_step=MTM_defined,
             name_alg='MTM')
  #### Examples with Rejection Free IIT    
    
    ex2_RF_IIT(theta=theta_selected,
               h=hsq_log,
               name_alg='IIT-RF-1',
               initial_K=.01) #tiny constany
    ############################    
    ex2_RF_IIT(theta=theta_selected,
               h=hsq_log,
               name_alg='IIT-RF-2',
               initial_K=1) #small constant  
    
    ############################
    ex2_RF_IIT(theta=theta_selected,
               h=hsq_log,
               name_alg='IIT-RF-4',
               initial_K=10*exp(theta_selected))#big bounding K
    
    ############################
    ex2_RF_IIT(theta=theta_selected,
               h=hsq_log,
               name_alg='IIT-RF-5',
               initial_K=300*exp(theta_selected))#huge bounding K    
    
  }
