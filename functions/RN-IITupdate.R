# Function to do 1 iteration of RN-IIT (Random Neighborhood IIT)
# We consider uniform proposal distribution over the neighbors
# Neighbor set includes all states that can be reached by swaping 1 coordinate
# But at each step we select m of those neighbors
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state) 
# h: balancing function
# m: Size of the Random Neighborhood
# Output: New state choosen proportionally and previous state's weight
RN_IITupdate <- function(X,pi,h,m){
  neighbors <- sample(1:p,m, replace = F) #Choose m neighbors
  probs <- numeric(m) #Vector to store weights
  pi_current <- pi(X)
  for(c in 1:m){
    Xnew <- X
    Xnew[neighbors[c]] <- 1-X[neighbors[c]] #Swap coordinate c
    probs[c] <- h(pi(Xnew)/pi_current)
  }
  
  U = runif(m)
  D = -log(U)/probs #Selecting a state proportional to the Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  Xnew[neighbors[index]] <- 1-X[neighbors[index]] #Choose neighbor state
  weight <- 1/mean(probs) #Compute weight using mean since we're using uniform dist.
  return(list(Xnew,weight))
}

### Second function 
### This one works using log probabilities and
### adjusting the balancing function accordingly
### In case probabilities are exp {something} and the balancing function 
### can be easily adapted
RN_IITupdate_log <- function(X,logpi,logh){
  neighbors <- sample(1:p,m, replace = F) #Choose m neighbors
  logprobs <- numeric(m) #Vector to store weights
  logpi_current <- logpi(X)
  for(c in 1:m){
    Xnew <- X
    Xnew[neighbors[c]] <- 1-X[neighbors[c]] #Swap coordinate c
    logprobs[c] <- logh(logpi(Xnew)-logpi_current)
  }
  
  U = runif(m)
  D = log(-log(U))-logprobs #Selecting a state proportional to the log Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  Xnew[neighbors[index]] <- 1-X[neighbors[index]] #Choose neighbor state
  weight <- 1/mean(exp(logprobs)) #Compute weight using mean since we're using uniform dist.
  return(list(Xnew,weight))
}
