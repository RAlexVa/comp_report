source(file.path(getwd(),'functions','balancing_functions.R')) #Balancing functions
source(file.path(getwd(),'functions','VT-IITupdate.R')) #VT-IIT functions
source(file.path(getwd(),'functions','IITupdate.R')) #Functions for IIT update
#Simulation study on a multimodal example
#rm(list=ls())






VT_IIT <- function(num_sim=50,
                   max_iter=5*10^4,
                   p=200,
                   n=100,
                   h_func=logh_func,
                   h_temp=hsq_log,
                   method,
                   J,
                   delta){
  simulation_name <- paste0('method_',method,'_J_',J,'_delta_',delta,'.csv')
  ##### Defining model #####
  L <- matrix(rnorm(n*p),nrow=n,ncol=p)
  
  beta <- runif(1,4,6)*sqrt(log(p)/n)
  Y <- (L[,1]+L[,2]+L[,3])*beta + rnorm(n,mean=0,sd=0.5)
  # Introducing multi modality
  L[,4] <- L[,2]-L[,3]+rnorm(n,mean=0,sd=0.1)
  L[,5] <- L[,1]+L[,2]+L[,3]+L[,6]+L[,7]+rnorm(n,mean=0,sd=0.1)
  L[,8] <- L[,6]-L[,7]+rnorm(n,mean=0,sd=0.1)
  
  ##### Function PI #####
  
  pi.distribution <- function(X){
    #Need to have Y and data defined outside the function
    selected <- which(X==1)
    subdata <- as.data.frame(cbind(Y,L[,selected]))
    
    if(sum(X)==0){
      res <- lm(Y~., data=subdata)
    }else{
      res <- lm(Y~.-1, data=subdata)
    }
    return(logLik(res)[1])
    #return(list(list(logLik(res)[1],res$coefficients))
  }
  # Testing function PI
  # pi.distribution(X=c(1,1,1,0,0,0,0,0,rep(0,p-8)))
  # pi.distribution(X=c(1,1,0,1,0,0,0,0,rep(0,p-8)))
  # pi.distribution(X=c(1,0,1,1,0,0,0,0,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,1,1,0,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,1,0,1,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,0,1,1,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,0,0,0,rep(0,p-8)))
# Vector of inverse temperatures
  vec_temp <-  1/(1+(delta*((1:J)-1))) #inverse-temperature vector
  
  ##### Create matrix to store data #####
  trajectory <- matrix(NA,nrow=num_sim*max_iter,ncol=13)
  #Add to the trajectory
  #column 1: simulation number
  #column 2: First 8 columns of state visited
  #Column 3: sum of other columns
  #column 4: current temperature
  #column 5: assigned weight
  
  ##### Run simulation #####
  
  for(sim in 1:num_sim){
    #Initialize the chain
    X <- rep(0,p)
    #Initialize the temperature
    curr_temp <- 1 #Start at temperature 1
    #Initialize other parameters
    logpsi <- rep(0,J)
    phi <- rep(1,length(vec_temp))
    for(i in 1:max_iter){
      if(i%%1000==1){print(paste0('VT-IIT (',method,') J=',J,' D=',delta,' sim:',sim,' iter:',i))}
      temp_neigh <- (1:J)[-curr_temp]
      #start_t <- Sys.time()
      iter <- VT_IITupdate_log(X,
                               curr_temp,
                               logpi=pi.distribution,
                               logh_func=h_func,
                               logh_temp=h_temp,
                               p,
                               vec_temp,
                               logpsi,
                               temp_neigh)
      #Sys.time()-start_t
      #Register stuff
      row_index <- (max_iter)*(sim-1) + i
      type <- iter[[5]]
      
      trajectory[row_index,1] <- sim
      trajectory[row_index,2:9] <- X[1:8]
      trajectory[row_index,10] <- sum(X[-(1:8)])
      trajectory[row_index,11] <- curr_temp
      trajectory[row_index,12] <-type
      trajectory[row_index,13] <-iter[[2]]
      #Update psi    
      logpsi <- update_logpsi(logpsi,curr_temp,i,J+1)
      
      #Update chain
      if(type=='temp'){
        curr_temp <- iter[[3]]
      }
      if(type=='space'){
        X <- iter[[1]]
      }
    }
  }
  write.csv(trajectory,file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}

just_IIT <- function(num_sim=100,
                   max_iter=5*10^4,
                   p=200,
                   n=100,
                   h_func,
                   method,
                   J,
                   delta){
  simulation_name <- paste0('method_',method,'_J_',J,'_delta_',delta,'.csv')
  ##### Defining model #####
  L <- matrix(rnorm(n*p),nrow=n,ncol=p)
  
  beta <- runif(1,4,6)*sqrt(log(p)/n)
  Y <- (L[,1]+L[,2]+L[,3])*beta + rnorm(n,mean=0,sd=0.5)
  # Introducing multi modality
  L[,4] <- L[,2]-L[,3]+rnorm(n,mean=0,sd=0.1)
  L[,5] <- L[,1]+L[,2]+L[,3]+L[,6]+L[,7]+rnorm(n,mean=0,sd=0.1)
  L[,8] <- L[,6]-L[,7]+rnorm(n,mean=0,sd=0.1)
  
  ##### Function PI #####
  
  pi.distribution <- function(X){
    #Need to have Y and data defined outside the function
    selected <- which(X==1)
    subdata <- as.data.frame(cbind(Y,L[,selected]))
    
    if(sum(X)==0){
      res <- lm(Y~., data=subdata)
    }else{
      res <- lm(Y~.-1, data=subdata)
    }
    return(logLik(res)[1])
    #return(list(list(logLik(res)[1],res$coefficients))
  }
  # Testing function PI
  # pi.distribution(X=c(1,1,1,0,0,0,0,0,rep(0,p-8)))
  # pi.distribution(X=c(1,1,0,1,0,0,0,0,rep(0,p-8)))
  # pi.distribution(X=c(1,0,1,1,0,0,0,0,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,1,1,0,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,1,0,1,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,0,1,1,rep(0,p-8)))
  # pi.distribution(X=c(0,0,0,0,1,0,0,0,rep(0,p-8)))
  
  
  ##### Create matrix to store data #####
  trajectory <- matrix(NA,nrow=num_sim*max_iter,ncol=11)
  #Add to the trajectory
  #column 1: simulation number
  #column 2: First 8 columns of state visited
  #column 3: assigned weight
  
  ##### Run simulation #####
  
  for(sim in 1:num_sim){
    #Initialize the chain
    X <- rep(0,p)
    #Initialize the temperature
    curr_temp <- 1 #Start at temperature 1
    for(i in 1:max_iter){
      if(i%%1000==1){print(paste0('VT-IIT (',method,') J=',J,' D=',delta,' sim:',sim,' iter:',i))}
      temp_neigh <- (1:J)[-curr_temp]
      #start_t <- Sys.time()
      iter <- IITupdate_log(X,logpi=pi.distribution,logh=h_func,p)
      #Sys.time()-start_t
      #Register stuff
      row_index <- (max_iter)*(sim-1) + i

      trajectory[row_index,1] <- sim
      trajectory[row_index,2:9] <- X[1:8]
      trajectory[row_index,10] <- sum(X[-(1:8)])
      trajectory[row_index,11] <-iter[[2]]

      X <- iter[[1]]
    }
  }
  write.csv(trajectory,file=file.path(getwd(),'results',simulation_name) ,row.names = F)
}

##### Running simulations #####
#Use the same seed for each simulation so all problems solve the same problem
# using some prompts 
