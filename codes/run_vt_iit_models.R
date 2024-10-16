#if(!require('Rcpp')){install.packages('Rcpp')}
library(Rcpp)
library(RcppArmadillo)
setwd('..')
Rcpp::sourceCpp("functions/cpp_functions.cpp")

#Define parameters for the model
seed_def <- 6055
n <- 100
p <- 200
iterations <- 50000
##### Choose chunk of simulations to run #####
writeLines('1 is simulations 1-20\n2 is simulations 21-40\n3 is simulations 41-60\n4 is simulations 61-80\n5 is simulations 81-100')
chunk_selected <- as.numeric(readline('Select chunk of simulations'))
if(chunk_selected<1 | chunk_selected>5){print('Error');chunk_selected <- -1;}
# chunk_selected <- 5
start_point <- 1+(chunk_selected-1)*20
end_point <- chunk_selected*20
start_point
end_point
#We should use the same seed to compare between models but each chunk should have a different seed.
#If we don't change the seed for each chunk the algorithm will follow the same trajectory for model 21 than for model 1
#although the likelihood would be different but still want to modify the random numbers used for sampling
set.seed(seed_def + chunk_selected)
##### Choose the temperature ladder to use #####
temp1.1 <- (1+((1:5)-1))^(-1) #J=4, Delta=1
temp1.2 <- (1+((1:10)-1))^(-1) #J=9, Delta=1
temp2.1 <- (1+((1:5)-1)*2)^(-1) #J=4, Delta=2
temp2.2 <- (1+((1:10)-1)*2)^(-1) #J=10, Delta=2

writeLines('1 is D=1,5 temperatures\n2 is D=1,10 temperatures\n3 is D=2,5 temperatures\n4 is D=2,10 temperatures')
t_selected <- as.numeric(readline('Select temperature ladder'))
if(t_selected==1){
  temp <- temp1.1  
}else if(t_selected==2){
  temp <- temp1.2
}else if(t_selected==3){
  temp <- temp2.1
}else if(t_selected==4){
  temp <- temp2.2
}else{print('Incorrect model selected')}

##### Choose method #####
#0 is iit, and there are other 3 methods depending on balancing functions for
writeLines('method 0 is IIT\nmethod 1 min bal. fun. for all\nmethod 2 half & half bal. fun.\nmethod 3 only max temp is different')
m_selected <- as.numeric(readline('Select method'))
if(m_selected==0){
  results <- Simulation_mod_IIT(n=n,p=p,startsim=start_point,endsim=end_point,numiter=iterations)  
}else if(m_selected==1){
  results <- Simulation_mod1(n=n,p=p,startsim=start_point,endsim=end_point,numiter=iterations,temp=temp,t=length(temp))  
}else if(m_selected==2){
  results <- Simulation_mod2(n=n,p=p,startsim=start_point,endsim=end_point,numiter=iterations,temp=temp,t=length(temp))
}else if(m_selected==3){
  results <- Simulation_mod3(n=n,p=p,startsim=start_point,endsim=end_point,numiter=iterations,temp=temp,t=length(temp))
}else{print('Incorrect model selected')}


write.table(results,paste0('results/','resultados_VT-IIT_modelo',m_selected,'_temp_',t_selected,'_seed_',seed_def,'+',chunk_selected,'sim',start_point,'_',end_point,'.csv'),row.names=F, col.names=F, sep=',')
