########### Testing the multimodality of the models ##############
#rm(list=ls())
library(Rcpp);library(RcppArmadillo);
Rcpp::sourceCpp("functions/cpp_testing_functions.cpp")
#Rcpp::sourceCpp("functions/cpp_functions.cpp")
mod1 <- c(1,1,1, rep(0,197))
mod2 <- c(1,1,0,1, rep(0,196))
mod3 <- c(1,0,1,1, rep(0,196))
mod4 <- c(0,0,0,0,1,1,1,rep(0,193))
mod5 <- c(0,0,0,0,1,1,0,1,rep(0,192))
mod6 <- c(0,0,0,0,1,0,1,1,rep(0,192))

check_modes(mod1)+1
check_modes(mod2)+1
check_modes(mod3)+1
check_modes(mod4)+1
check_modes(mod5)+1
check_modes(mod6)+1



########### Checking multi modality in all models ##############
modes <- list(mod1,mod2,mod3,mod4,mod5,mod6) #All modes in a vector
comparison <- array(-1,dim=c(100,6,200))
for(model_selected in 1:100){
  print(paste0('model ',model_selected))
  Y_res <- read_Y_cpp(paste0('models/resY',model_selected,'.csv'))
  X_model <- read_file_cpp(paste0('models/modelX',model_selected,'.csv'))
  for(m in 1:length(modes)){#loop for modes
    tocheck <- modes[[m]]
    a <- logLikelihood(X_model, Y_res,which(tocheck==1)-1)
    for(i in 1:length(tocheck)){
      newstate <- tocheck
      newstate[i] <- 1-newstate[i]
      if(logLikelihood(X_model, Y_res,which(newstate==1)-1)>a){
        # print(paste0('happened in ',model_selected,"_",m,"_",i))
        # print(paste0("likelihoods mode:",round(a,3)," other:",round(logLikelihood(X_model, Y_res,which(newstate==1)-1),2)))
        comparison[model_selected,m,i] <- 1
      }
      
    }
  }
}
#which(comparison==-1, arr.ind=T)
which(comparison==1, arr.ind=T) #147 models with issues
comparison[3,3,2]
Y_res <- read_Y_cpp(paste0('models/resY',3,'.csv'))
X_model <- read_file_cpp(paste0('models/modelX',3,'.csv'))
logLikelihood_m(X_model,Y_res,which(mod3==1)-1)
logLikelihood_m(X_model,Y_res,which(c(1,1,1,1)==1)-1)

########### Checking individual modes ##############

selected <- 100
Y_res <- read_Y_cpp(paste0('models/resY',selected,'.csv'))
X_model <- read_file_cpp(paste0('models/modelX',selected,'.csv'))

logLikelihood_m(X_model, Y_res,which(mod1==1)-1)
logLikelihood_m(X_model, Y_res,which(mod2==1)-1)
logLikelihood_m(X_model, Y_res,which(mod3==1)-1)
logLikelihood_m(X_model, Y_res,which(mod4==1)-1)
logLikelihood_m(X_model, Y_res,which(mod5==1)-1)
logLikelihood_m(X_model, Y_res,which(mod6==1)-1)

tocheck <- mod1
tocheck <- c(1,1,1,rep(0,123),1,rep(0,53))
tocheck <- c(1,1,1,rep(0,123),1,rep(0,23),1,rep(0,49))
a <- logLikelihood(X_model, Y_res,which(tocheck==1)-1)
compare <- c()
for(i in 1:length(tocheck)){
  newstate <- tocheck
  newstate[i] <- 1-newstate[i]
  compare[i] <- logLikelihood(X_model, Y_res,which(newstate==1)-1)>a
}



tocheck <- mod3
a <- logLikelihood(X_model, Y_res,which(tocheck==1)-1)
compare <- c()
for(i in 1:length(tocheck)){
  newstate <- tocheck
  newstate[i] <- 1-newstate[i]
  compare[i] <- logLikelihood(X_model, Y_res,which(newstate==1)-1)>a
}

comparison[1,2,3]
Y_res <- read_Y_cpp(paste0('models/resY',1,'.csv'))
X_model <- read_file_cpp(paste0('models/modelX',1,'.csv'))
logLikelihood_m(X_model,Y_res,which(mod2==1)-1)
logLikelihood_m(X_model,Y_res,which(c(1,1,1,1)==1)-1)


tocheck <- mod3
a <- logLikelihood(X_model, Y_res,which(tocheck==1)-1)
compare <- c()
for(i in 1:length(tocheck)){
  newstate <- tocheck
  newstate[i] <- 1-newstate[i]
  compare[i] <- logLikelihood(X_model, Y_res,which(newstate==1)-1)>a
}