# Function to do 1 iteration of MH-IIT
# We consider uniform proposal distribution over the neighbors
# Neighbors are all states that can be reached by swaping 1 coordinate
# Input:
# pi: Function to return value of unnormalized target distribution
# X: Vector of size p (Current state) 
# h: balancing function
# rho: proportion to use IIT or MH
# Output: New state choosen proportionally and previous state's weight
source(file.path(getwd(),'functions','IITupdate.R')) #Functions for IIT update
source(file.path(getwd(),'functions','MHupdate.R')) #Functions for MH update

MH_IITupdate <- function(X,pi,h,p,rho){
  if(runif(1)<=rho){
    return(IITupdate(X,pi,h,p))
  }else{
    return(MHupdate(X,pi,h,p))
  }
}

### Second function 
### This one works using log probabilities and
### adjusting the balancing function accordingly
### In case probabilities are exp {something} and the balancing function 
### can be easily adapted
MH_IITupdate_log <- function(X,logpi,logh,p,rho){
  if(runif(1)<=rho){
    return(IITupdate_log(X,logpi,logh,p))
  }else{
    return(MHupdate_log(X,logpi,logh,p))
  }
}
