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
temp <- temp2.1

test_1 <- Simulation_mod3(n,p,1, 10,100, temp,length(temp))
test_2 <- Simulation_mod_IIT(n,p,1, 10,100)

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
#rm(list=ls())
Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")
#Rcpp::sourceCpp("functions/cpp_functions.cpp")
check_modes(c(1,1,1,rep(0,197)))+1
check_modes(c(1,1,0,1,rep(0,196)))+1
check_modes(c(1,0,1,1,rep(0,196)))+1
check_modes(c(0,0,0,0,1,1,1,rep(0,193)))+1
check_modes(c(0,0,0,0,1,1,0,1,rep(0,192)))+1
check_modes(c(0,0,0,0,1,0,1,1,rep(0,192)))+1

check_modes(c(1,1,1,1,rep(0,196)))+1
########### Checking multimodality ##############
Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")
mod1 <- c(1,1,1, rep(0,197))
mod2 <- c(1,1,0,1, rep(0,196))
mod3 <- c(1,0,1,1, rep(0,196))
mod4 <- c(0,0,0,0,1,1,1,rep(0,193))
mod5 <- c(0,0,0,0,1,1,0,1,rep(0,192))
mod6 <- c(0,0,0,0,1,0,1,1,rep(0,192))

selected <- 18
Y_res <- read_Y_cpp(paste0('models/resY',selected,'.csv'))
X_model <- read_file_cpp(paste0('models/modelX',selected,'.csv'))

logLikelihood_m(X_model, Y_res,which(mod1==1)-1)
logLikelihood_m(X_model, Y_res,which(mod2==1)-1)
logLikelihood_m(X_model, Y_res,which(mod3==1)-1)
logLikelihood_m(X_model, Y_res,which(mod4==1)-1)
logLikelihood_m(X_model, Y_res,which(mod5==1)-1)
logLikelihood_m(X_model, Y_res,which(mod6==1)-1)

tocheck <- mod6
(a <- logLikelihood_m(X_model, Y_res,which(tocheck==1)-1))

compare <- c()
for(i in 1:length(mod1)){
  newstate <- tocheck
  newstate[i] <- 1-newstate[i]
  lll <- logLikelihood(X_model, Y_res,which(newstate==1)-1)
  print(lll)
  compare[i] <- lll>a
}
sum(compare)
#Checking this likelihood
mod5 <- c(0,0,0,0,1,1,0,1,rep(0,192))
logLikelihood(X_model, Y_res,which(mod5==1)-1)
mod5_comp <- c(0,0,0,0,1,1,0,1,rep(0,152),1,rep(0,39))
logLikelihood(X_model, Y_res,which(mod5_comp==1)-1)
compare <- c()
for(i in 1:length(mod5)){
  newstate <- mod5
  newstate[i] <- 1-newstate[i]
  compare[i] <- logLikelihood(X_model, Y_res,which(newstate==1)-1)
}

logLikelihood(X_model, Y_res,which(c(1,1,1,1,1,1,1,1,1,1,1,1,1)==1)-1)

logLikelihood(X_model, Y_res,which(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)==1)-1)

logLikelihood(X_model, Y_res,which(c(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1)==1)-1)
#Computing likelihood in R
variables <- c(1:3,9:20)
variables <- c(1:3)
variables <- which(c(1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1)==1)
data <- as.data.frame(cbind(Y_res,X_model))
data_s <- data[,c(1,variables+1)] #Include column 1 since it's Y, shift indexes to consider column Y
mod <- lm(V1 ~.-1, data=data_s)
logLik(mod)
#Computing likelihood using Rcpp functions
logLikelihood(X_model,Y_res,variables-1) #Shift the indexes

###################### Checking if the normals are properly created ##################
normals <- gen_normals(1000)
hist(normals$n1) #Small variance
hist(normals$n2) #Large variance
hist(normals$n3) #Mean=5
hist(normals$n4) #standard normal





########### Testing if the functions are working properly reading the file ##############

Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")
Rcpp::sourceCpp("functions/cpp_functions.cpp")
llik <- read_file_inside_function(100)
Simulation_mod1(n=n,p=p,startsim=1,endsim=100,numiter=1,temp=temp,t=length(temp))
Simulation_mod2(n=n,p=p,startsim=1,endsim=100,numiter=1,temp=temp,t=length(temp))
Simulation_mod3(n=n,p=p,startsim=1,endsim=100,numiter=1,temp=temp,t=length(temp))
Simulation_mod_IIT(n=n,p=p,startsim=1,endsim=100,numiter=1)


########### Testing IIT function que está escribiendo arameo ##############
Rcpp::sourceCpp("functions/cpp_functions.cpp")
n <- 100
p <- 200
test1 <- Simulation_mod_IIT(n=n,p=p,startsim = 1, endsim=20,numiter=1000)


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

######  Testing functions for parallel tempering #############
#rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions.cpp")
Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")

test_loglik(c(1,1,1),matX,resY)
logLikelihood(matX,resY,c(0,1,2))
test_loglik(c(0,0,0,0,0),matX,resY)
logL_0(resY)
# RF_PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap, vec temp, SEXP method_input)

test <- RF_PT_IIT_sim(p=200,startsim=1, endsim=1, numiter=50, iterswap=10, temp=1/1:5, method='M1')


matX <- readmodelX(1)
resY <- readY(1)
cur <- logL_0(resY)
Z_factor(rep(0,200),'sq',modelX=matX,resY=resY, temperature=1)
Z_factor(rep(0,200),'min',modelX=matX,resY=resY, temperature=1)
Z_factor(c(1,1,1,rep(0,197)),'sq',modelX=matX,resY=resY, temperature=1) #Small Z_h for a mode
Z_factor(c(1,1,1,rep(0,197)),'min',modelX=matX,resY=resY, temperature=1) #Small Z_h for a mode
# Z_factor(vec X, String chosen_bf, mat modelX, vec resY, double temperature)
logLikelihood(matX,resY,c(0))-cur
logLikelihood(matX,resY,c(1))-cur
logLikelihood(matX,resY,c(2))-cur
logLikelihood(matX,resY,c(3))-cur
logLikelihood(matX,resY,c(4))-cur
logLikelihood(matX,resY,c(5))-cur
logLikelihood(matX,resY,c(6))-cur

RF_IIT_sim(2,100,200)

test_coord(c(1,0,0,0,0,1,0,0,1))
test_coord(c(0,0,0,0,0))

temp1.1 <- (1+((1:5)-1))^(-1)
temp1.2 <- (1+((1:10)-1))^(-1)
temp2.1 <- (1+((1:5)-1)*2)^(-1)
temp2.2 <- (1+((1:10)-1)*2)^(-1)
temp <- temp1.1

check <- RF_PT_IIT_sim(p=200,startsim=1,endsim=1,numiter=100,iterswap=200,temp=temp,method="M1")

ccc <- RF_update(rep(0,20), "sq",matX,resY)

########### Testing the index process ##############
#rm(list=ls())
Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")
#Rcpp::sourceCpp("functions/cpp_functions.cpp")
random_binom(5,3,1,2)
set.seed(123)
random_binom(7,5,1,2)


#index_pro(int p, vec temp, int numswap, double swap_prob )
index_pro(10,c(1,1/2,1/3,1/4,1/5),20,0.7)
