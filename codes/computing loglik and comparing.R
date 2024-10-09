library(Rcpp)

Rcpp::sourceCpp("functions/cpp_functions.cpp")
#Rcpp::sourceCpp("codes/loglikelihood.cpp")

########################################
#Below we compare the computation of likelihood between R and Cpp to see it matches
res <- random_model(n=100,p=200)
#Computing likelihood in R
variables <- c(1:3,9:20)
data <- as.data.frame(cbind(res$Y,res$X))
data_s <- data[,c(1,variables+1)] #Include column 1 since it's Y, shift indexes to consider column Y
mod <- lm(V1 ~.-1, data=data_s)
logLik(mod)
#Computing likelihood using Rcpp functions
logLikelihood(res$X,res$Y,variables-1) #Shift the indexes

###############################################
#Now computing likelihood without features
variables <- c()
data_s <- as.data.frame(cbind(res$Y,res$X[,variables+1]))
mod <- lm(V1 ~., data=data_s)
logLik(mod)

###### Testing how to redefine a function by fixing some parameters
logL_0(res$Y)



cppFunction("double logLikelihood(uvec pos) {
   bool result = (num % 2 == 1);
   return result;
}")






set.seed(234)
lm_r <- function(n,p){
  res <- random_model(n,p)
  data <- as.data.frame(cbind(res$Y,res$X))
  data_s <- data[,c(1:4,10:20)]
  mod <- lm(V1 ~.-1, data=data_s)
  lll <- logLik(mod)
  return(lll[1])
}
set.seed(234)
lm_cpp <- function(n,p){
  res <- random_model(n,p)
  data <- as.data.frame(cbind(res$Y,res$X))
  data_s <- data[,c(1:4,10:20)]
  mod <- fastLm(as.matrix(data_s[,-1]),as.numeric(unlist(data_s[1])))
  #Computing log-likelihood
  #https://www.statlect.com/fundamentals-of-statistics/linear-regression-maximum-likelihood
  
  return(logLikelihood_num(mod$residuals))
}

set.seed(234)
lm_r(100,200)
set.seed(234)
lm_cpp(100,200)
microbenchmark(lm_r(100,200),lm_cpp(100,200))


