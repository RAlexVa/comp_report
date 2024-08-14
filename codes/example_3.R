setwd('..')
#rm(list=ls())

### Source functions
source(file.path(getwd(),'functions','IITupdate.R')) #Functions for IIT update
source(file.path(getwd(),'functions','MHupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','MH-IITupdate.R')) #Functions for MH update
source(file.path(getwd(),'functions','RN-IITupdate.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','MTM.R')) #Functions for RN-IIT update
source(file.path(getwd(),'functions','balancing_functions.R')) #Balancing functions
source(file.path(getwd(),'functions','IIT-RFupdate.R')) #Functions for IIT update
example3 <- function(num_sim=50,
                     max_iter=500*1000,
                     p=200,
                     p_1,
                     theta,
                     threshold=0.5,
                     h,
                     update_step,
                     name_alg){
  simulation_name <- paste0('ex3_',name_alg,'_t',theta,'_p',p,'_p1',p_1,'_sim',num_sim,'.csv')
  ### Call and define functions to use for the simulation
  X_mode1 <- c(1,0,rep(1,p_1-1),rep(0,p-p_1-1)) #Global mode for example 2
  X_mode2 <- c(0,rep(1,p_1),rep(0,p-p_1-1))
  #
  pi.distribution <- function(X){ #Function to return the  log-probabilities
    d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
    d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
    
    maxd <- max(d1,d2)
    mind <- min(d1,d2)
    #We use log(A+B) = log(A) + log(1 + B/A) approx= log(A) + B/A (If B/A is small)
    expres <- -theta*mind + log(1+exp(-theta*(maxd-mind)))
    return(expres)
  }
  
  # pi.d2 <- function(X){
  #   d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
  #   d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
  #   expres <- exp(-theta*d1) + exp(-theta*d2)
  #   return(expres)
  # }
  
  Fm <- function(X){
    d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
    d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
    return(c(d1,d2))
  }
  # True pi_F
  pi_F_true <- matrix(0,nrow=p+1,ncol=p+1)
  for(n in 0:(p-2)){
    coord <- n+1
    #First coordinates are 10
    pi_F_true[coord,coord+2] <- choose(p-2,n) * (exp(-theta*n) + exp(-theta*(n+2)))
    #First coordinates are 01
    pi_F_true[coord+2,coord] <- choose(p-2,n) * (exp(-theta*n) + exp(-theta*(n+2)))
    #First coordinates are 00 and 11 
    pi_F_true[coord+1,coord+1] <- choose(p-2,n)*4*exp(-theta*(n+1))
  }
  pi_F_true <- pi_F_true/(2*(1+exp(-theta))^p)
  
  # To measure convergence
  dist_pi <- function(est){ #Function to return the TVD between the true pi and the estimated pi (pushforward measures)
    if(!identical(dim(est),dim(pi_F_true))){
      print('Size of matrix is not the same')
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
    pi_F_est <- matrix(0,nrow=p+1,ncol=p+1)
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
      Feval <- Fm(X) +c(1,1) #Function F evaluated
      pi_F_est[Feval[1],Feval[2]] <- pi_F_est[Feval[1],Feval[2]] + W #Assign estimated weight
      
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
ex3_RN_IIT <- function(num_sim=50,
                     max_iter=500*1000,
                     p=200,
                     p_1,
                     theta,
                     threshold=0.5,
                     h,
                     update_step=RN_IITupdate_log,
                     m=40,
                     name_alg){
  simulation_name <- paste0('ex3_',name_alg,'_t',theta,'_p',p,'_p1',p_1,'_sim',num_sim,'.csv')
  ### Call and define functions to use for the simulation
  X_mode1 <- c(1,0,rep(1,p_1-1),rep(0,p-p_1-1)) #Global mode for example 2
  X_mode2 <- c(0,rep(1,p_1),rep(0,p-p_1-1))
  #
  pi.distribution <- function(X){ #Function to return the  log-probabilities
    d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
    d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
    
    maxd <- max(d1,d2)
    mind <- min(d1,d2)
    #We use log(A+B) = log(A) + log(1 + B/A) approx= log(A) + B/A (If B/A is small)
    expres <- -theta*mind + log(1+exp(-theta*(maxd-mind)))
    return(expres)
  }
  
  # pi.d2 <- function(X){
  #   d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
  #   d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
  #   expres <- exp(-theta*d1) + exp(-theta*d2)
  #   return(expres)
  # }
  
  Fm <- function(X){
    d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
    d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
    return(c(d1,d2))
  }
  # True pi_F
  pi_F_true <- matrix(0,nrow=p+1,ncol=p+1)
  for(n in 0:(p-2)){
    coord <- n+1
    #First coordinates are 10
    pi_F_true[coord,coord+2] <- choose(p-2,n) * (exp(-theta*n) + exp(-theta*(n+2)))
    #First coordinates are 01
    pi_F_true[coord+2,coord] <- choose(p-2,n) * (exp(-theta*n) + exp(-theta*(n+2)))
    #First coordinates are 00 and 11 
    pi_F_true[coord+1,coord+1] <- choose(p-2,n)*4*exp(-theta*(n+1))
  }
  pi_F_true <- pi_F_true/(2*(1+exp(-theta))^p)
  
  # To measure convergence
  dist_pi <- function(est){ #Function to return the TVD between the true pi and the estimated pi (pushforward measures)
    if(!identical(dim(est),dim(pi_F_true))){
      print('Size of matrix is not the same')
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
    pi_F_est <- matrix(0,nrow=p+1,ncol=p+1)
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
      #Adding a minimum that is above 0 to avoid issues with very low numbers
      Feval <- Fm(X) +c(1,1) #Function F evaluated
      pi_F_est[Feval[1],Feval[2]] <- pi_F_est[Feval[1],Feval[2]] + W #Assign estimated weight
      
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
ex3_RF_IIT <- function(num_sim=50,
                     max_iter=500*1000,
                     p=200,
                     p_1,
                     theta,
                     threshold=0.5,
                     h,
                     update_step=IIT_RFupdate_log,
                     name_alg,
                     initial_K){
  simulation_name <- paste0('ex3_',name_alg,'_t',theta,'_p',p,'_p1',p_1,'_sim',num_sim,'_K',round(initial_K,2),'.csv')
  ### Call and define functions to use for the simulation
  X_mode1 <- c(1,0,rep(1,p_1-1),rep(0,p-p_1-1)) #Global mode for example 2
  X_mode2 <- c(0,rep(1,p_1),rep(0,p-p_1-1))
  #
  pi.distribution <- function(X){ #Function to return the  log-probabilities
    d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
    d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
    
    maxd <- max(d1,d2)
    mind <- min(d1,d2)
    #We use log(A+B) = log(A) + log(1 + B/A) approx= log(A) + B/A (If B/A is small)
    expres <- -theta*mind + log(1+exp(-theta*(maxd-mind)))
    return(expres)
  }
  
  # pi.d2 <- function(X){
  #   d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
  #   d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
  #   expres <- exp(-theta*d1) + exp(-theta*d2)
  #   return(expres)
  # }
  
  Fm <- function(X){
    d1 <- sum(abs(X-X_mode1)) #Distance to mode 1
    d2 <- sum(abs(X-X_mode2)) #Distance to mode 2
    return(c(d1,d2))
  }
  # True pi_F
  pi_F_true <- matrix(0,nrow=p+1,ncol=p+1)
  for(n in 0:(p-2)){
    coord <- n+1
    #First coordinates are 10
    pi_F_true[coord,coord+2] <- choose(p-2,n) * (exp(-theta*n) + exp(-theta*(n+2)))
    #First coordinates are 01
    pi_F_true[coord+2,coord] <- choose(p-2,n) * (exp(-theta*n) + exp(-theta*(n+2)))
    #First coordinates are 00 and 11 
    pi_F_true[coord+1,coord+1] <- choose(p-2,n)*4*exp(-theta*(n+1))
  }
  pi_F_true <- pi_F_true/(2*(1+exp(-theta))^p)
  
  # To measure convergence
  dist_pi <- function(est){ #Function to return the TVD between the true pi and the estimated pi (pushforward measures)
    if(!identical(dim(est),dim(pi_F_true))){
      print('Size of matrix is not the same')
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
    pi_F_est <- matrix(0,nrow=p+1,ncol=p+1)
    early_finish <- FALSE #To check if convergence was achieved before max. iterations
    bounding_K <- initial_K #Specific for RF-IIT
    #Initialize
    X <- rep(0,p)
    #Counting the number of times we used the PI function
    count_PIs <- 0
    #Counting the number of times we updated K
    count_K <- 0 #Specific for RF-IIT
    ### Within that for loop need a loop for the steps (considering max number of iterations)
    for(step in 1:max_iter){
      if(step %% 1000 ==0){print(paste(name_alg,'theta',theta,',Iteration:',step,',Simulation:',i))}
      iter <- update_step(X,pi.distribution,h,p,bounding_K)
      W <- max(iter[[2]],.Machine$double.eps) #Estimated weight of previous state, bounding it from below
      count_PIs <- count_PIs+iter[[3]]
      bounding_K <- iter[[4]]#Specific for RF-IIT
      count_K <- count_K+iter[[5]]#Specific for RF-IIT
      #Adding a minimum that is above 0 to avoid issues with very low numbers
      Feval <- Fm(X) +c(1,1) #Function F evaluated
      pi_F_est[Feval[1],Feval[2]] <- pi_F_est[Feval[1],Feval[2]] + W #Assign estimated weight
      
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
    iteration_K <- c(iteration_K,count_K)#Specific for RF-IIT
    final_K <- c(final_K,bounding_K)#Specific for RF-IIT
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
    
    example3(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             update_step=IITupdate_log,
             name_alg='IIT')
    
    example3(theta=theta_selected,p_1=p1_selected,
             h=hmin_log,
             update_step=MHupdate_log,
             name_alg='MH')
    
    MH_IIT_defined <- function(X,pi,h,p){
      return(MH_IITupdate_log(X,pi,h,p,rho=0.025))
    }
    example3(theta=theta_selected,p_1=p1_selected,
             h=hmin_log,
             update_step=MH_IIT_defined,
             name_alg='MH-IIT')
    
    ex3_RN_IIT(theta=theta_selected,p_1=p1_selected,
               h=hmin_log,
               name_alg='RN-IIT')
    
    MTM_defined <- function(X,pi,h,p){
      return(MTMupdate_log(X,pi,h,p,m=40))
    }
    example3(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             update_step=MTM_defined,
             name_alg='MTM')
#### Examples with Rejection Free IIT    

    ex3_RF_IIT(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-1',
             initial_K=.01) #tiny constany
    ############################    
    ex3_RF_IIT(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-2',
             initial_K=1) #small constant  

    ############################
    ex3_RF_IIT(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-4',
             initial_K=10*exp(theta_selected))#big bounding K
    
    ############################
    ex3_RF_IIT(theta=theta_selected,p_1=p1_selected,
             h=hsq_log,
             name_alg='IIT-RF-5',
             initial_K=300*exp(theta_selected))#huge bounding K    
  }
}




