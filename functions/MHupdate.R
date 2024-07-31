# Function to do 1 iteration of uninformed Metropolis Hastings
# We consider uniform proposal distribution over the neighbors
# Neighbors are all states that can be reached by swaping 1 coordinate
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state) 
# Output: New state 
MHupdate <- function(X,pi,h){
  c <- sample(1:p,size=1) #Choose a coordinate at random
  Xnew[c] <- 1-X[c] #Swap coordinate c
  prob <- min(1,pi(Xnew)/pi(X)) #Compute probability of accepting
  U = runif(1)
  if(U>=prob){ #If rejected
    Xnew <- X #We stay in the current state
  }
  weight <- 1#The weight is always 1 
  return(list(Xnew,weight))
}

### Second function 
### This one works using log probabilities 
### In case probabilities are exp {something}
MHupdate_log <- function(X,logpi,logh){
  c <- sample(1:p,size=1) #Choose a coordinate at random
  Xnew[c] <- 1-X[c] #Swap coordinate c
  logprob <- min(0,logpi(Xnew)-logpi(X)) #Compute probability of accepting
  U = runif(1)
  if(log(U)>=logprob){ #If rejected
    Xnew <- X #We stay in the current state
  }
  weight <- 1 #The weight is always 1 
  return(list(Xnew,weight))
}
