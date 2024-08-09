### Source functions
#rm(list=ls())
setwd('..')
source(file.path(getwd(),'functions','IIT-RFupdate.R')) #Functions for IIT update
source(file.path(getwd(),'functions','balancing_functions.R')) #Balancing functions

example1 <- function(num_sim=50,
                     max_iter=500*1000,
                     p=500,
                     p_1,
                     theta,
                     threshold=0.1,
                     h,
                     update_step=IIT_RFupdate_log,
                     name_alg,
                     initial_K){
  simulation_name <- paste0('ex1_',name_alg,'_t',theta,'_p',p,'_p1_',p_1,'_sim',num_sim,'K',round(initial_K,2),'.csv')
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
  iteration_K <- c()
  final_K <- c()
  for(i in 1:num_sim){
    print(paste(name_alg,'Simulation:',i,'theta',theta))
    ### To check thresholds
    pi_F_est <- numeric(p+1)
    early_finish <- FALSE #To check if convergence was achieved before max. iterations
    bounding_K <- initial_K 
    #Initialize
    X <- rep(0,p)
    #Counting the number of times we used the PI function
    count_PIs <- 0
    #Counting the number of times we updated K
    count_K <- 0
    ### Within that for loop need a loop for the steps (considering max number of iterations)
    for(step in 1:max_iter){
      if(step %% 1000 ==0){print(paste(name_alg,'theta',theta,',Iteration:',step,',Simulation:',i))}
      iter <- update_step(X,pi.distribution,h,p,bounding_K)
      W <- max(iter[[2]],.Machine$double.eps) #Estimated weight of previous state, bounding it from below
      count_PIs <- count_PIs+iter[[3]]
      bounding_K <- iter[[4]]
      count_K <- count_K+iter[[5]]
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
    iteration_K <- c(iteration_K,count_K)
    final_K <- c(final_K,bounding_K)
    #If after max_iter iterations it didn't converge
    if(!early_finish){iter_conv <- c(iter_conv,max_iter)} #register max_iter  
  }
  #iter_conv
  write.csv(cbind(iter_conv,calls_for_pi,last_F_value,iteration_K,final_K),file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}
# Parameters
thetas <- c(6,7,8,9) 
p1s <- c(20,50)
set.seed(348) #Define seed
for(i in 1:length(thetas)){
  for(j in 1:length(p1s)){
    theta_selected <- thetas[i]
    p1_selected <- p1s[j]
############################    
    example1(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-1',
             initial_K=.01) #timy constany
############################    
    example1(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-2',
             initial_K=1) #minimum constant  
############################    
    example1(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-3',
             initial_K=exp(theta_selected/2))#Optimal bounding K
    
############################
    example1(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-4',
             initial_K=10*exp(theta_selected))#big bounding K
    
    ############################
    example1(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-5',
             initial_K=300*exp(theta_selected))#huge bounding K
  }
}
