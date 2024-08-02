# Function to do 1 iteration of MTM (Multiple Try MH with some balancing function for the weights)
# We consider uniform proposal distribution over the neighbors
# Neighbor set includes all states that can be reached by swaping 1 coordinate
# But at each step we select m of those neighbors
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state) 
# h: balancing function (The weight is h applied to the ratio of probabilities)
# m: Number of proposals to consider
# Output: New state choosen proportionally and previous state's weight
MTMupdate <- function(X,pi,h,p,m){
  neighbors <- sample(1:p,m, replace = F) #Choose m DISTINCT neighbors
  probs_going <- numeric(m) #Vector to store weights
  pi_current <- pi(X)
  for(c in 1:m){
    Xnew <- X
    Xnew[neighbors[c]] <- 1-X[neighbors[c]] #Swap coordinate c
    probs_going[c] <- h(pi(Xnew)/pi_current)
  }
  probs_going <- probs_going/sum(probs_going) #Normalize the probabilities
  Y <- X
  choose_neigh <- sample(1:m,size=1,prob=probs_going) #Choose which neighbor
  Y[neighbors[choose_neigh]] <- 1-Y[neighbors[choose_neigh]] #Swap coordinate to choose a neighbor Y
  
  # X is also a neighbor of Y so we restrict the selection to not include it
  # It's manually added later on to the possible neighbors
  Yneighbors <- sample((1:p)[-neighbors[choose_neigh]],m,replace=F) #Choose m DISTINCT neighbors from selected Y
  Yneighbors[choose_neigh] <- neighbors[choose_neigh] #Add X as a neighbor
  probs_returning <- numeric(m) #Vector to store weights
  piY_current <- pi(Y)
  for(c in 1:m){
    Ynew <- Y
    Ynew[Yneighbors[c]] <- 1-Y[Yneighbors[c]] #Swap coordinate c
    probs_returning[c] <- h(pi(Ynew)/piY_current)
  }
  probs_returning <- probs_returning/sum(probs_returning)
  
  Ap <- (piY_current/pi_current)*(probs_returning[choose_neigh]/probs_going[choose_neigh])#computing accetance probability
  Ap <- min(1,Ap)
  if(runif(1)<Ap){
    Xnew <- Y
  }else{
    Xnew <- X
  }
  weight <- 1 #Weight is always 1
  return(list(Xnew,weight))
}

### Second function 
### This one works using log probabilities and
### adjusting the balancing function accordingly
### In case probabilities are exp {something} and the balancing function 
### can be easily adapted
MTMupdate_log <- function(X,logpi,logh,p,m){
  neighbors <- sample(1:p,m, replace = F) #Choose m neighbors
  logprobs_going <- numeric(m) #Vector to store weights
  logpi_current <- logpi(X)
  for(c in 1:m){
    Xnew <- X
    Xnew[neighbors[c]] <- 1-X[neighbors[c]] #Swap coordinate c
    logprobs_going[c] <- logh(logpi(Xnew)-logpi_current)
  }
  probs_going <- exp(logprobs_going) #Normalize the probabilities
  probs_going <- probs_going/sum(probs_going) #Normalize the probabilities
  Y <- X
  choose_neigh <- sample(1:m,size=1,prob=probs_going) #Choose which neighbor
  Y[neighbors[choose_neigh]] <- 1-Y[neighbors[choose_neigh]] #Swap coordinate to choose a neighbor Y
  
  # X is also a neighbor of Y so we restrict the selection to not include it
  # It's manually added later on to the possible neighbors
  Yneighbors <- sample((1:p)[-neighbors[choose_neigh]],m,replace=F) #Choose m DISTINCT neighbors from selected Y
  Yneighbors[choose_neigh] <- neighbors[choose_neigh] #Add X as a neighbor
  logprobs_returning <- numeric(m) #Vector to store weights
  logpiY_current <- logpi(Y)
  for(c in 1:m){
    Ynew <- Y
    Ynew[Yneighbors[c]] <- 1-Y[Yneighbors[c]] #Swap coordinate c
    logprobs_returning[c] <- logh(logpi(Ynew)-logpiY_current)
  }
  probs_returning <- exp(logprobs_returning)
  probs_returning <- probs_returning/sum(probs_returning)
  
  Ap <- exp(logpiY_current-logpi_current)*(probs_returning[choose_neigh]/probs_going[choose_neigh])#computing accetance probability
  Ap <- min(1,Ap)
  if(runif(1)<Ap){
    Xnew <- Y
  }else{
    Xnew <- X
  }
  weight <- 1 #Weight is always 1
  return(list(Xnew,weight))
}
