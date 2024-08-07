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
  return(list(Xnew,weight,m+m-1))
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
  # print('new iteration')
  # print(logpi_current)
  # print(sum(X))
  # print(X)
  # print(neighbors)
  # print(X[neighbors])
  # print(logprobs_going)
  probs_going <- exp(logprobs_going) #Apply exp to the probabilities
  probs_going[probs_going==0] <- .Machine$double.eps; #Fix extremely low probabilities
  coord <- which(probs_going==Inf)
  if(!identical(integer(0),coord)){#Fix extremely high probabilities
    probs_going[coord] <- 1;
    probs_going[-coord] <- .Machine$double.eps;
  }
  #print(probs_going)
  probs_going <- probs_going/sum(probs_going) #Normalize the probabilities
  #print(probs_going)
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
  probs_returning[probs_returning==0] <- .Machine$double.eps; #Fix extremely low probabilities
  coord <- which(probs_returning==Inf)
  if(!identical(integer(0),coord)){#Fix extremely high probabilities
    probs_returning[coord] <- 1;
    probs_returning[-coord] <- .Machine$double.eps;
  }
  probs_returning <- probs_returning/sum(probs_returning)
  # print('new iteration')
  # print(probs_returning)
  # print(logpiY_current)
  # print(logpi_current)
  # print(choose_neigh)
  # print(probs_returning[choose_neigh])
  # print(probs_going[choose_neigh])
  Ap <- exp(logpiY_current-logpi_current)*(probs_returning[choose_neigh]/probs_going[choose_neigh])#computing acceptance probability
  Ap <- min(1,Ap)
  if(runif(1)<Ap){
    Xnew <- Y
  }else{
    Xnew <- X
  }
  weight <- 1 #Weight is always 1
  return(list(Xnew,weight,m+m-1))
}
