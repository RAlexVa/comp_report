#rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
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
#########################################################
##### Testing functiosn that read model files instead of creating them each time ######
Rcpp::sourceCpp("functions/cpp_functions.cpp")
n <- 100
p <- 200
set.seed(134)
temp1.1 <- (1+((1:5)-1))^(-1)
temp1.2 <- (1+((1:10)-1))^(-1)
temp2.1 <- (1+((1:5)-1)*2)^(-1)
temp2.2 <- (1+((1:10)-1)*2)^(-1)
temp <- temp1.1
test_1 <- Simulation_mod1(n=n,p=p,startsim=1,endsim=100,numiter=30,temp=temp,t=length(temp))



######################################################
#Checking results
library(data.table)
res1 <- fread('results/resultados_VT-IIT_modelo1_temp_1_seed_123.csv')
res2 <- fread('results/resultados_VT-IIT_modelo2_temp_1_seed_123.csv')

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

####Checking latest function
Rcpp::sourceCpp("functions/cpp_functions.cpp")
n <- 100
p <- 200
set.seed(134)
temp1.1 <- (1+((1:5)-1))^(-1)
temp1.2 <- (1+((1:10)-1))^(-1)
temp2.1 <- (1+((1:5)-1)*2)^(-1)
temp2.2 <- (1+((1:10)-1)*2)^(-1)
temp <- temp1.1
test_1 <- Simulation_mod1(n=n,p=p,numsim=2,numiter=1001,temp=temp,t=length(temp))

test_2 <- Simulation_mod2(n=n,p=p,numsim=3,numiter=50,temp=temp,t=length(temp))

test_3 <- Simulation_mod3(n=n,p=p,numsim=5,numiter=50,temp=temp,t=length(temp))

test_4 <- Simulation_mod_IIT(n=n,p=p,numsim=5,numiter=50)


n <- 100
p <- 200
set.seed(134)
before <- Sys.time()
test_1 <- Simulation_mod1_full(n=n,p=p,numsim=1,numiter=1000,temp=temp,t=length(temp))
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



Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")
#check how to read files from C++
# Create many models
set.seed(345)
for(i in 1:20){
  mod <- random_model(n=100,p=200)
  write.table(mod$X,file=paste0('models/modelX',i,'.csv'),sep=',',row.names = F,col.names=F);
  write.table(mod$Y,paste0('models/resY',i,'.csv'),sep=',',row.names = F,col.names=F);
}

mmm <- read_file_cpp('models/modelX1.csv')

selected <- 2
Y_res <- read_Y_cpp(paste0('models/resY',selected,'.csv'))
X_model <- read_file_cpp(paste0('models/modelX',selected,'.csv'))
logLikelihood(X_model,Y_res,c(0,1,3))

Y_res <- read_file_inside_function(5)


read_file_inside_function(3)


check_ones(c(1,0,0,0,0,1,0,1,1,0,0,1,0,1))
check_modes(c(1,0,1,0))

check_modes(c(1,0,1,0,0,0,0,1,1,1,1))

check_modes(c(1,1,1,0,0,0,0,0,0,1,1,1,1))

########### Testing the specific modes ##############
check_modes(c(1,1,1,rep(0,197)))+1
check_modes(c(1,1,0,1,rep(0,196)))+1
check_modes(c(1,0,1,1,rep(0,196)))+1
check_modes(c(0,0,0,0,1,1,1,rep(0,193)))+1
check_modes(c(0,0,0,0,1,1,0,1,rep(0,192)))+1
check_modes(c(0,0,0,0,1,0,1,1,rep(0,192)))+1





########### Testing if the functions are working properly reading the file ##############

Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")
Rcpp::sourceCpp("functions/cpp_functions.cpp")
llik <- read_file_inside_function(100)
Simulation_mod1(n=n,p=p,startsim=1,endsim=100,numiter=1,temp=temp,t=length(temp))
Simulation_mod2(n=n,p=p,startsim=1,endsim=100,numiter=1,temp=temp,t=length(temp))
Simulation_mod3(n=n,p=p,startsim=1,endsim=100,numiter=1,temp=temp,t=length(temp))
Simulation_mod_IIT(n=n,p=p,startsim=1,endsim=100,numiter=1)





#######################################
sample_prop(c(4,2,3,4))


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
