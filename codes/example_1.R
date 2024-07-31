### Source functions
source(file.path(getwd(),'functions','IITupdate.R')) #Functions for IIT update
source(file.path(getwd(),'functions','balancing_functions.R')) #Balancing functions


example1 <- function(num_sim, max_iter, p,p_1,theta,threshold,h,update_step){
  simulation_name <- paste0('ex1_t',theta,'_p',p,'_p1',p_1,'.csv')
  ### Call and define functions to use for the simulation
  X_mode <- c(rep(1,p_1),rep(0,p-p_1)) #Global mode for example 1
  # pi <- function(X){ #Function to return the unnormalized target pi
  #   dif <- -sum(abs(X-X_mode))
  #   return(exp(theta*dif))
  # }
  pi <- function(X){ #Function to return the log probabilities
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
  for(i in 1:num_sim){
    print(paste('Simulation:',i))
    ### To check thresholds
    pi_F_est <- numeric(p+1)
    early_finish <- FALSE #To check if convergence was achieved before max. iterations
    #Initialize
    X <- rep(0,p)
    ### Within that for loop need a loop for the steps (considering max number of iterations)
    for(step in 1:max_iter){
      if(step %% 1000 ==0){print(paste('Iteration: ',step,', Simulation:',i))}
      #  iter <- update_step(X,pi,h)
      iter <- update_step(X,pi,h)
      W <- iter[[2]] #Estimated weight of previous state
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
    #If after max_iter iterations it didn't converge
    if(!early_finish){iter_conv <- c(iter_conv,max_iter)} #register max_iter  
  }
  #iter_conv
  write.csv(iter_conv,file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}

set.seed(348) #Define seed
example1(num_sim=50,
         max_iter=500*1000,
         p=500,
         p_1=50,
         theta=6,
         threshold=0.1,
         h=hmin_log,
         update_step=IITupdate_log)
