library(Rcpp)

Rcpp::sourceCpp("codes/random_normal.cpp")
Rcpp::sourceCpp("codes/loglikelihood.cpp")

res <- random_model(n=100,p=200)
data <- as.data.frame(cbind(res$Y,res$X))
#data_s <- data
data_s <- data[,c(1:4,10:20)]

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


