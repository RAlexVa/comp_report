# Function to do 1 iteration of IIT with Rejection Free
# We consider uniform proposal distribution over the neighbors
# Neighbors are all states that can be reached by swaping 1 coordinate
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state) 
# h: balancing function
# Output: New state choosen proportionally, previous state's weight, number of calls for pi, updated bounding constant
IIT_RFupdate <- function(X,pi,h,p, bounding_K){
  K <- bounding_K #Use global variable K
  probs <- numeric(p) #Vector to store relative importance 
  pi_current <- pi(X)
  for(c in 1:p){
    Xnew <- X
    Xnew[c] <- 1-X[c] #Swap coordinate c
    probs[c] <- h(pi(Xnew)/pi_current)
  }
  newK <- max(K,max(probs),1/min(probs)) #Update the bounding constant K
  newK <- min(newK,1.797693e+308) #Bound the max so we dont get Inf
  U = runif(p)
  D = -log(U)/probs #Selecting a state proportional to the Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  Xnew[index] <- 1-X[index] #Choose neighbor state
  param <- max(mean(probs)/newK,.Machine$double.xmin)
  weight <- 1 + rgeom(1,param) #Compute weight using a geometric distribution
  while(is.na(weight)){weight <- 1 + rgeom(1,param)} #In case of NA, run again
  return(list(Xnew,weight,p,newK,(K!=newK)*1))
}

### Second function 
### This one works using log probabilities and
### adjusting the balancing function accordingly
### In case probabilities are exp {something} and the balancing function 
### can be easily adapted
IIT_RFupdate_log <- function(X,logpi,logh,p,bounding_K){
  K <- bounding_K #Use global variable K
  logprobs <- numeric(p) #Vector to store weights
  logpi_current <- logpi(X)
  for(c in 1:p){
    Xnew <- X
    Xnew[c] <- 1-X[c] #Swap coordinate c
    logprobs[c] <- logh(logpi(Xnew)-logpi_current)
  }
  newK <- max(K,exp(max(logprobs)),exp(-min(logprobs))) #Update the bounding constant K
  newK <- min(newK,.Machine$double.xmax) #Bound the max so we dont get Inf
  U = runif(p)
  D = log(-log(U))-logprobs #Selecting a state proportional to the log Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  Xnew[index] <- 1-X[index] #Choose neighbor state
  param <- max(mean(exp(logprobs))/newK,.Machine$double.xmin)
  weight <- 1 + rgeom(1,param) #Compute weight using a geometric distribution
  while(is.na(weight)){weight <- 1 + rgeom(1,param)} #In case of NA, run again
  return(list(Xnew,weight,p,newK,(K!=newK)*1))
}
