#rm(list=ls())
library(tidyverse)
library(latex2exp)
files <- dir(file.path(getwd(),'results_moba'))

# Remove files of only 1 simulation
# sim1 <- files[grep('sim1.csv',files)]
# 
# for(f in sim1){
#   file.remove(file.path(getwd(),'results_moba',f))
# }

#ex1 <- files[grep('ex1_',files)]
# rf_iit <- tibble(ex,alg,theta,p_1,ini_K,iter_conv,calls_for_pi,last_F_value,iteration_K,final_K)
# data_ex1 <- tibble(ex,alg,theta,p_1,iter_conv,calls_for_pi,last_F_value)
# data_ex2 <- tibble(ex,alg,theta,iter_conv,calls_for_pi,last_F_value)

ex1_rf_iit <- tibble()
data_ex1 <- tibble() 
data_ex2 <- tibble()

for(i in 1:length(files)){
  temp_file <- read_csv(file.path(getwd(),'results_moba',files[i]), col_types = cols())
  parameters <- unlist(strsplit(gsub('[sim,t,p]','',gsub('.csv','',files[i])),split="[_,K]"))
  if(parameters[2]=='IIT-RF' & parameters[1]=='ex1'){
    par_matrix <- matrix(parameters[c(1,2,3,6,8)],nrow=50,ncol=5,byrow=T)
    colnames(par_matrix) <- c('ex','alg','theta','p_1','ini_K')
    ex1_rf_iit <- rbind(ex1_rf_iit,cbind(par_matrix,temp_file))
  }
  if(parameters[1]=='ex1' & parameters[2]!='IIT-RF'){
    par_matrix <- matrix(parameters[c(1,2,3,6)],nrow=50,ncol=4,byrow=T)
    colnames(par_matrix) <- c('ex','alg','theta','p_1')
    data_ex1 <- rbind(data_ex1,cbind(par_matrix,temp_file))
  }
  if(parameters[1]=='ex2' & parameters[2]!='IIT-RF' ){
    par_matrix <- matrix(parameters[c(1,2,3)],nrow=50,ncol=3,byrow=T)
    colnames(par_matrix) <- c('ex','alg','theta')
    data_ex2 <- rbind(data_ex2,cbind(par_matrix,temp_file))
  }
}
#Fix dataset for ex1 with RF-IIT
fix <- unique(ex1_rf_iit[,c('p_1','theta','ini_K')])
fix <- cbind(fix,'k'=rep(1:5,nrow(fix)/5))
colnames(fix)[length(fix)] <- 'K'
ex1_rf_iit <- left_join(ex1_rf_iit,fix,by=c('p_1','theta','ini_K'))
ex1_rf_iit$alg <- paste0(ex1_rf_iit$alg,'-',ex1_rf_iit$K)
#Add it to the ex1 dataset
##### Plot for example 1 #####

data_ex1 <- rbind(data_ex1,ex1_rf_iit |> select(colnames(data_ex1)))

data_ex1 |> 
  filter(p_1==50) |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 1')

data_ex1 |> 
  filter(p_1==20) |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 1')

##### Plot for example 2 #####

data_ex2 |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  filter(alg!='RN-IIT') |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))

