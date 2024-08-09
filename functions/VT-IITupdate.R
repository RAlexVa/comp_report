# Function to do 1 iteration of VT-IIT (Varying temperature IIT)
# We consider uniform proposal distribution over the neighbors
# Neighbor set includes all states that can be reached by swaping 1 coordinate
# We consider uniform proposal distribution over the temperatures 
# We can modify the neighbor set for temperatures
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state)
# vec_temp: Vector of size J+1 (contains all inverse temperatures in decreasing order)
# curr_temp: index (Current temperature)
# temp_neigh: vector indicating the index of temperatures to consider as neighbors
#             by default is all temperatures except the current one.
# h_func: list with balancing functions to use when moving through states
# h_temp: balancing function to use when moving through temperatures
# inv_temp: vector with inverse temperatures to use
# phi: Vector to use in the algorithm
# psi: Vector to use in the algorithm
# Output: New state choosen proportionally and previous state's weight
VT_IITupdate <- function(X,curr_temp,pi,h_func,h_temp,p,vec_temp,psi,phi=rep(1,length(vec_temp)),temp_neigh=(1:length(vec_temp))[-curr_temp],n_0=100,s_0=100){
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
    temp_weight <- h_temp((pi_current^(vec_temp[t]-temp_new)) * (psi[t]/psi[curr_temp]))
    temprobs <- c(temprobs,temp_weight)
  }
  J <- length(temp_neigh)
  joint_probs <- c(probs/p,temprobs/J) #Multiply by the respective Q and get all in 1 vector
  U = runif(length(joint_probs))
  D = -log(U)/joint_probs #Selecting a state proportional to the Prob. vector
  index = which.min(D) #Select index to choose neighbor
  
  Xnew <- X
  if(index<=p){#If a space neighbor is selected
    Xnew[index] <- 1-X[index] #Choose neighbor state
    newtemp <- curr_temp #don't change temperature
  }if(index>p){#If a temperature neighbor is selected
    Xnew <- X #Don't change state
    newtemp <- index-p
  }
  

#Compute weight
  inv_w <- pi_current^(1-temp_now)*phi(curr_temp)/(psi(curr_temp)*sum(joint_probs))
  weight <- 1/inv_w #Compute weight using mean since we're using uniform dist.
  
#Update PSI (According to the paper)
  iteration <- i #This comes from outside
  
#decrease current psi
  temp <- psi[curr_temp]*exp(-s_0/(n_0 + iteration))
#Increase other psi
  psi <-  psi[curr_temp]*exp(s_0/((J-1)*(n_0 + iteration))) #J counts the number of temperatures, including 0
  psi[curr_temp] <- temp
  
  return(list(Xnew,weight,newtemp,psi,p+J))
}

### Second function 
### This one works using log probabilities and
### adjusting the balancing function accordingly
### In case probabilities are exp {something} and the balancing function 
### can be easily adapted

