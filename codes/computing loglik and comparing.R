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



######################################################
#This was in the cpp_functions file before
# library(Rcpp)
Rcpp::sourceCpp("functions/cpp_functions.cpp")
#rm(list=ls())
n <- 10
p <- 15
set.seed(134)
#res <- random_model(n,p)

temp <- c(1,1.2,2.3,3.5,3.9,4)
logpsi <- c(1,2,3,4,5,6)
logpsi2 <- rep(0,6)
length(temp)==length(logpsi)
set.seed(134)

n <- 100
p <- 200
set.seed(134)
before <- Sys.time()
test_1 <- Simulation_mod1(n=n,p=p,numsim=1,numiter=1000,temp=temp,t=length(temp))
Sys.time()-before
dim(test_1$states)
dim(test_1$temps)
dim(test_1$logpsi)


#Modes that 
visited <- test_1$states

#Checks
table(check_temp)
which(check_temp!=0)[1:3] #Where temperature changed from 0
which(check_temp>1)[1:3] #Where temperature reached 2
check_temp[29:40]
rowSums(visited[29:40,])
identical(visited[29,],visited[30,])
identical(visited[30,],visited[31,]) #Ok since we changed temperature from 0 to 1
identical(visited[31,],visited[32,]) #Ok since we changed temperature from 1 to 0
identical(visited[32,],visited[33,]) #Ok since we didn't change temperature

X <- c(0,1,0,rep(0,12))

curr_temp <- 3
set.seed(134)
(check <- VT_IIT_update_c1(X=X,
                           temp=temp,
                           curr_temp=curr_temp,
                           modelX=res$X,
                           resY=res$Y,
                           n=length(X),
                           t=length(temp),
                           logpsi=logpsi)) #number of temperatures



test_vec(10)

ma <- matrix(0,nrow=p,ncol=p)
X_input <- rep(0,p);
#mat test_mat(mat M, vec X, vec t, int curr_temp, int p)
test_mat(ma,X_input,p)

# update_logpsi(c(1,2,3,4,5),2,10,5)

# library(dplyr)
# set.seed(134)
# temperature <- temp[curr_temp+1]
# coord <- which(X==1)+1
# t <- length(temp)
# logpi_r <- c(p+t)
# data <- as.data.frame(cbind(res$Y,res$X))
# logpi_current <- logLik(lm(V1 ~.-1, data=as.data.frame(data[,c(1,coord)])))[1]
# for(i in 1:p){
#   Xnew <- X;
#   Xnew[i] <- 1-X[i];
#   variables <- which(Xnew==1);
#   if(length(variables)==0){
#     data_s <- data |>  select('V1');
#     mod <- lm(V1 ~., data=data_s)
#   }else{
#     data_s <- data[,c(1,variables+1)] #Include column 1 since it's Y, shift indexes to consider column Y
#     mod <- lm(V1 ~.-1, data=data_s)
#   }
# logpi_r[i] <- exp(min((logLik(mod)[1] - logpi_current)*temperature,0))/length(X)
# }
# for(i in 1:t){
#   if(temp[i]==temperature){logpi_r[p+i] <- 0}else{
#   logpi_r[p+i] <- exp(min(logpi_current*(temp[i]-temperature),0))/(length(temp)-1)
# }}
# logpi_r
# u <- runif(p+t)
# which.min(-log(u)/logpi_r)


### Compare
# fun_dummy_1 <- function(X=X,
#                         temp=temp,
#                         curr_temp=curr_temp,
#                         modelX=res$X,
#                         resY=res$Y,
#                         n=length(X),
#                         t=length(temp)){
#   return(VT_IIT_update_c1(X=X,
#                           temp=temp,
#                           curr_temp=curr_temp,
#                           modelX=res$X,
#                           resY=res$Y,
#                           n=length(X),
#                           t=length(temp)))
# }
# 
# fun_dummy_2 <- function(X=X,
#                         temp=temp,
#                         curr_temp=curr_temp,
#                         modelX=res$X,
#                         resY=res$Y,
#                         n=length(X),
#                         t=length(temp)){
#   temperature <- temp[curr_temp+1]
#   coord <- which(X==1)+1
#   t <- length(temp)
#   logpi_r <- c(p+t)
#   data <- as.data.frame(cbind(res$Y,res$X))
#   logpi_current <- logLik(lm(V1 ~.-1, data=as.data.frame(data[,c(1,coord)])))[1]
#   for(i in 1:p){
#     Xnew <- X;
#     Xnew[i] <- 1-X[i];
#     variables <- which(Xnew==1);
#     if(length(variables)==0){
#       data_s <- data |>  select('V1');
#       mod <- lm(V1 ~., data=data_s)
#     }else{
#       data_s <- data[,c(1,variables+1)] #Include column 1 since it's Y, shift indexes to consider column Y
#       mod <- lm(V1 ~.-1, data=data_s)
#     }
#     logpi_r[i] <- exp(min((logLik(mod)[1] - logpi_current)*temperature,0))/length(X)
#   }
#   for(i in 1:t){
#     if(temp[i]==temperature){logpi_r[p+i] <- 0}else{
#       logpi_r[p+i] <- exp(min(logpi_current*(temp[i]-temperature),0))/(length(temp)-1)
#     }}
#   logpi_r
#   
#   u <- runif(p+t)
#   
#   index_n <- which.min(-log(u)/logpi_r)
#   if(index_n<=p){X[index_n] <- 1-X[index_n]}else{
#     temperature = temp[index_n-p]
#   }
#   
#   return(list(temperature,logpi_r,X,))
# }
# 
# microbenchmark(fun_dummy_1, fun_dummy_2, times=1000)
