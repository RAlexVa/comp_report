#install.packages('tidyverse')
#install.packages('Rcpp')

n <- 100
p <- 200
set.seed(234)
L <- matrix(rnorm(n*p),nrow=n,ncol=p)

beta <- runif(1,4,6)*sqrt(log(p)/n)
Y <- (L[,1]+L[,2]+L[,3])*beta + rnorm(n,mean=0,sd=0.5)
# Introducing multi modality
L[,4] <- L[,2]-L[,3]+rnorm(n,mean=0,sd=0.1)
L[,5] <- L[,1]+L[,2]+L[,3]+L[,6]+L[,7]+rnorm(n,mean=0,sd=0.1)
L[,8] <- L[,6]-L[,7]+rnorm(n,mean=0,sd=0.1)

##### Function PI #####

pi.distribution <- function(X){
  #Need to have Y and data defined outside the function
  selected <- which(X==1)
  subdata <- as.data.frame(cbind(Y,L[,selected]))
  
  if(sum(X)==0){
    res <- lm(Y~., data=subdata)
  }else{
    res <- lm(Y~.-1, data=subdata)
  }
  return(logLik(res)[1])
  #return(list(list(logLik(res)[1],res$coefficients))
}

X <- rep(0,p)

pi.distribution(X)

n <- 100
p <- 200
set.seed(234)
L <- matrix(rnorm(n*p),nrow=n,ncol=p)

beta <- runif(1,4,6)*sqrt(log(p)/n)
Y <- (L[,1]+L[,2]+L[,3])*beta + rnorm(n,mean=0,sd=0.5)
# Introducing multi modality
L[,4] <- L[,2]-L[,3]+rnorm(n,mean=0,sd=0.1)
L[,5] <- L[,1]+L[,2]+L[,3]+L[,6]+L[,7]+rnorm(n,mean=0,sd=0.1)
L[,8] <- L[,6]-L[,7]+rnorm(n,mean=0,sd=0.1)

# Manually extracting coefficients and log-likelihood information
X <- c(rep(1,3),rep(p-3))
selected <- which(X==1)
subdata <- as.data.frame(cbind(Y,L[,selected]))
mod <- lm(Y~.-1, data=subdata)

M <- as.matrix(subdata[,-1])
coef <- solve(t(M)%*%M)%*%t(M)%*%Y
res <- Y-M%*%coef #residuals
sq_res <- as.numeric(t(res)%*%res) #Sum of squared residuals
# sig <- sqrt(t(res)%*%res/97) #sigma of residuals
# s2 <- as.numeric(t(res)%*%res/97)
# summary(mod)$sigma
# t(res)%*%res*(1/(dim(M)[1])) #

s2 <- as.numeric(t(res)%*%res/100)
logl <- (n*log(2*pi) + n*log(s2) + sq_res/s2)/(-2)

(n*log(2*pi) + n*log(as.numeric(t(res)%*%res)) -n*log(n) + n)/(-2)
