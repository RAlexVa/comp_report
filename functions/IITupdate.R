# Function to do 1 iteration of IIT (Informed Importance Tempering)
# We consider uniform proposal distribution over the neighbors
# Neighbors are all states that can be reached by swaping 1 coordinate
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state) 
# h: balancing function
# Output: New state choosen proportionally and previous state's weight
IITupdate <- function(X,pi,h){
  probs <- numeric(p) #Vector to store weights
  pi_current <- pi(X)
  for(c in 1:p){
    Xnew <- X
    Xnew[c] <- 1-X[c] #Swap coordinate c
    probs[c] <- h(pi(Xnew)/pi_current)
  }
  
  U = runif(p)
  D = -log(U)/probs #Selecting a state proportional to the Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  Xnew[index] <- 1-X[index] #Choose neighbor state
  weight <- 1/mean(probs) #Compute weight using mean since we're using uniform dist.
  return(list(Xnew,weight))
}

### Second function 
### This one works using log probabilities and
### adjusting the balancing function accordingly
### In case probabilities are exp {something} and the balancing function 
### can be easily adapted
IITupdate_log <- function(X,logpi,logh){
  logprobs <- numeric(p) #Vector to store weights
  logpi_current <- logpi(X)
  for(c in 1:p){
    Xnew <- X
    Xnew[c] <- 1-X[c] #Swap coordinate c
    logprobs[c] <- logh(logpi(Xnew)-logpi_current)
  }
  
  U = runif(p)
  D = log(-log(U))-logprobs #Selecting a state proportional to the log Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  Xnew[index] <- 1-X[index] #Choose neighbor state
  weight <- 1/mean(exp(logprobs)) #Compute weight using mean since we're using uniform dist.
  return(list(Xnew,weight))
}
