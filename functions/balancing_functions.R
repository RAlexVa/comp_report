### Balancing functions

hmin <- function(r){return(min(1,r))}
hsq <- function(r){return(sqrt(r))}
hmax <- function(r){return(max(1,r))}
hc <- function(r,c=1){
  first <- min(1,r*exp(-c))
  second <- min(r,exp(-c))
  return(max(first,second))
}


### Version for log probabilities
hmin_log <- function(r){return(min(0,r))}
hsq_log <- function(r){return(r/2)}
hmax_log <- function(r){return(max(0,r))}
hc_log <- function(r,c=1){
  first <- min(0,r-c)
  second <- min(r,-c)
  return(max(first,second))
}