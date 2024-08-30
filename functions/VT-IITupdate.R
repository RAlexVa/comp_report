# Function to do 1 iteration of VT-IIT (Varying temperature IIT)
# We consider uniform proposal distribution over the neighbors
# Neighbor set includes all states that can be reached by swaping 1 coordinate
# We consider uniform proposal distribution over the temperatures 
# We can modify the neighbor set for temperatures
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state)
# vec_temp: Vector of size J (contains all inverse temperatures in decreasing order)
# curr_temp: index (Current temperature)
# temp_neigh: vector indicating the index of temperatures to consider as neighbors
#             by default is all temperatures except the current one.
# h_func: list with balancing functions to use when moving through states
# h_temp: balancing function to use when moving through temperatures
# inv_temp: vector with inverse temperatures to use
# phi: In this case we just use 1 but it should be a Vector with same size as temperatures to use in the algorithm
# psi: Vector to use in the algorithm
# Output: New state choosen proportionally and previous state's weight
VT_IITupdate <- function(X,curr_temp,pi,h_func,h_temp,p,vec_temp,psi,temp_neigh,phi=1){
temp_now <- vec_temp[curr_temp] #Temperature now
h <- h_func[[curr_temp]] #Balancing function to use depending on the temperature
#Computing weight
#Checking all spacial neighbors
  probs <- numeric(p) #Vector to store weights
  pi_current <- pi(X)
  for(c in 1:p){
    Xnew <- X
    Xnew[c] <- 1-X[c] #Swap coordinate c
    probs[c] <- h((pi(Xnew)/pi_current)^(temp_now))
  }
  
#Checking all temperature neighbors
  temprobs <- c()
  for(t in temp_neigh){
    temp_weight <- h_temp((pi_current^(vec_temp[t]-temp_now)) * (psi[t]/psi[curr_temp]))
    temprobs <- c(temprobs,temp_weight)
  }
  J <- length(temp_neigh)
  joint_probs <- c(probs/p,temprobs/J) #Multiply by the respective Q and get all in 1 vector
  U = runif(length(joint_probs))
  D = -log(U)/joint_probs #Selecting a state proportional to the Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  if(index<=p){#If a space neighbor is selected
    Xnew[index] <- 1-X[index] #Update state
    newtemp <- curr_temp #Don't change temperature
    jump <- 'space' #Indicate movement in state
  }
  if(index>p){#If a temperature neighbor is selected
    Xnew <- X #Don't change state
    newtemp <- temp_neigh[index-p] #Update temperature index
    jump <- 'temp' #Indicate movement in temp
  }
  

#Compute weight
  #weight <- pi_current^(1-temp_now)*phi[curr_temp]/(psi[curr_temp]*sum(joint_probs))
  weight <- pi_current^(1-temp_now)/(psi[curr_temp]*sum(joint_probs))
  
  return(list(Xnew,weight,newtemp,p+J,jump))
}

update_psi <- function(psi,curr_temp,iteration,J,n_0=100,s_0=100){
  #Update PSI (According to the paper)
  #decrease current psi
  temp <- psi[curr_temp]*exp(-s_0/(n_0 + iteration))
  #Increase other psi
  psi <-  psi*exp(s_0/((J-1)*(n_0 + iteration))) #J counts the number of temperatures, including 0
  psi[curr_temp] <- temp
  return(psi)
}

### Second function 
### This one works using log probabilities and
### adjusting the balancing function accordingly
### In case probabilities are exp {something} and the balancing function 
### can be easily adapted
VT_IITupdate_log <- function(X,curr_temp,logpi,logh_func,logh_temp,p,vec_temp,logpsi,temp_neigh,phi=1){
  temp_now <- vec_temp[curr_temp] #Temperature now
  logh <- logh_func[[curr_temp]] #Balancing function to use depending on the temperature
  #Computing weight
  #Checking all spacial neighbors
  logprobs <- numeric(p) #Vector to store weights
  logpi_current <- logpi(X)
  for(c in 1:p){
    Xnew <- X
    Xnew[c] <- 1-X[c] #Swap coordinate c
    logprobs[c] <- logh((logpi(Xnew)-logpi_current)*(temp_now))
  }
  
  #Checking all temperature neighbors
  logtemprobs <- c()
  for(t in temp_neigh){
    #Using psi
    #temp_weight <- logh_temp((logpi_current*(vec_temp[t]-temp_now)) * (log(psi[t])-log(psi[curr_temp])))
    #Using logpsi
    logtemp_weight <- logh_temp((logpi_current*(vec_temp[t]-temp_now)) + (logpsi[t]-logpsi[curr_temp]))
    logtemprobs <- c(logtemprobs,logtemp_weight)
  }
  J <- length(temp_neigh)
  #Usando logprobs
  joint_logprobs <- c(logprobs-log(p),logtemprobs-log(J)) #Multiply by the respective Q and get all in 1 vector
  U = runif(length(joint_logprobs))
  D = log(-log(U))-joint_logprobs #Selecting a state proportional to the log Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  if(index<=p){#If a space neighbor is selected
    Xnew[index] <- 1-X[index] #Update state
    newtemp <- curr_temp #Don't change temperature
    jump <- 'space' #Indicate movement in state
  }
  if(index>p){#If a temperature neighbor is selected
    Xnew <- X #Don't change state
    newtemp <- temp_neigh[index-p] #Update temperature index
    jump <- 'temp' #Indicate movement in temp
  }
  
  
  #Compute weight
  #weight <- exp(logpi_current*(1-temp_now))*phi[curr_temp]/(psi[curr_temp]*sum(joint_probs))
  weight <- exp(logpi_current*(1-temp_now)-logpsi[curr_temp])/(sum(exp(joint_logprobs)))
  return(list(Xnew,weight,newtemp,p+J,jump))
}

update_logpsi <- function(logpsi,curr_temp,iteration,J,n_0=100,s_0=100){
  #Update PSI (According to the paper)
  #decrease current psi
  temp <- logpsi[curr_temp] - s_0/(n_0 + iteration)
  #Increase other psi
  logpsi <-  logpsi + s_0/((J-1)*(n_0 + iteration)) #J counts the number of temperatures, including 0
  logpsi[curr_temp] <- temp
  return(logpsi)
}
