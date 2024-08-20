#rm(list=ls())
library(tidyverse)
library(latex2exp)
library(Hmisc)
file_route <- file.path('results_moba','3_run_14AGO')
file_route <- file.path('results_moba')
files <- dir(file.path(getwd(),file_route))
files <- files[grepl('.csv',files)] #Get only CSV files, ignore folders

ex1_rf_iit <- tibble()
ex2_rf_iit <- tibble()
ex3_rf_iit <- tibble()
data_ex1 <- tibble() 
data_ex2 <- tibble()
data_ex3 <- tibble()
##### Read all files #####
for(i in 1:length(files)){
  temp_file <- read_csv(file.path(getwd(),file_route,files[i]), col_types = cols())
  parameters <- unlist(strsplit(gsub('[sim,t,p]','',gsub('.csv','',files[i])),split="[_,K]"))
  if(parameters[1]=='ex1' & grepl('^IIT-RF',parameters[2])){
    par_matrix <- matrix(parameters[c(1,2,3,6,8)],nrow=50,ncol=5,byrow=T)
    colnames(par_matrix) <- c('ex','alg','theta','p_1','ini_K')
    ex1_rf_iit <- rbind(ex1_rf_iit,cbind(par_matrix,temp_file))
  }
  if(parameters[1]=='ex1' & !grepl('^IIT-RF',parameters[2])){
    par_matrix <- matrix(parameters[c(1,2,3,6)],nrow=50,ncol=4,byrow=T)
    colnames(par_matrix) <- c('ex','alg','theta','p_1')
    data_ex1 <- rbind(data_ex1,cbind(par_matrix,temp_file))
  }
  if(parameters[1]=='ex2' & !grepl('^IIT-RF',parameters[2]) ){
    par_matrix <- matrix(parameters[c(1,2,3)],nrow=50,ncol=3,byrow=T)
    colnames(par_matrix) <- c('ex','alg','theta')
    data_ex2 <- rbind(data_ex2,cbind(par_matrix,temp_file))
  }
  if(parameters[1]=='ex3' & !grepl('^IIT-RF',parameters[2])){
    if(substr(parameters[5],1,1)==1){#Parameter p1 was not properly labeled
      parameters[5] <- substr(parameters[5],2,nchar(parameters[5]))
      par_matrix <- matrix(parameters[c(1,2,3,5)],nrow=50,ncol=4,byrow=T)
      colnames(par_matrix) <- c('ex','alg','theta','p_1')
      data_ex3 <- rbind(data_ex3,cbind(par_matrix,temp_file))
    }
  }
  if(parameters[1]=='ex3' & grepl('^IIT-RF',parameters[2])){
    if(substr(parameters[5],1,1)==1){#Parameter p1 was not properly labeled
      parameters[5] <- substr(parameters[5],2,nchar(parameters[5]))
      par_matrix <- matrix(parameters[c(1,2,3,5,8)],nrow=50,ncol=5,byrow=T)
      colnames(par_matrix) <- c('ex','alg','theta','p_1','ini_K')
      ex3_rf_iit <- rbind(ex3_rf_iit,cbind(par_matrix,temp_file))
    }
    
  }
}
##### Fix dataset for ex1 with RF-IIT #####

fix <- unique(ex1_rf_iit[,c('p_1','theta','ini_K')])
fix <- cbind(fix,'k'=rep(1:5,nrow(fix)/5))
colnames(fix)[length(fix)] <- 'K'
ex1_rf_iit <- left_join(ex1_rf_iit,fix,by=c('p_1','theta','ini_K'))
ex1_rf_iit$alg <- paste0(ex1_rf_iit$alg,'-',ex1_rf_iit$K)
#Add it to the ex1 dataset
data_ex1 <- rbind(data_ex1,ex1_rf_iit |> select(colnames(data_ex1)))


##### Fix dataset for ex3 with RF-IIT #####
#Add it to the ex1 dataset
data_ex3 <- rbind(data_ex3,ex3_rf_iit |> select(colnames(data_ex1)))
  

##### Plot for example 1 #####

ex1_img1 <- data_ex1 |> 
  filter(p_1==50) |> filter(alg %nin% paste0('IIT-RF-',1:5)) |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 1')
ex1_img1 #Show
#export
jpeg(file.path(getwd(),'plots',"ex1.jpg"), width = 700)
ex1_img1
dev.off() 


ex1_img2 <- data_ex1 |> 
  filter(p_1==50) |> filter(alg %in% c(paste0('IIT-RF-',1:5),'IIT')) |> 
  #filter(theta>8) |> 
  #mutate(post_calls=calls_for_pi) |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 1')
ex1_img2 #Show
#export
jpeg(file.path(getwd(),'plots',"ex1_2.jpg"), width = 700)
ex1_img2
dev.off() 


#export for slides
jpeg(file.path(getwd(),'plots',"ex1_2_67.jpg"), width = 700)
data_ex1 |> 
  filter(p_1==50) |> filter(alg %in% c(paste0('IIT-RF-',1:5),'IIT')) |> 
  filter(theta%in%c(6,7)) |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 1')
dev.off() 
jpeg(file.path(getwd(),'plots',"ex1_2_89.jpg"), width = 700)
data_ex1 |> 
  filter(p_1==50) |> filter(alg %in% c(paste0('IIT-RF-',1:5),'IIT')) |> 
  filter(theta%in%c(8,9)) |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 1')
dev.off() 



data_ex1 |> 
  filter(p_1==20) |>  filter(alg!='MTM', alg %nin% paste0('IIT-RF-',1:5)) |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 1')

##### Plot for example 2 #####

ex2 <- data_ex2 |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  filter(alg%nin%c('IIT')) |> 
  filter(alg%nin%c('MTM')) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 2')
ex2 #Show
#export
jpeg(file.path(getwd(),'plots',"ex2.jpg"), width = 700)
ex2
dev.off() 

##### Plot for example 3 #####
ex3 <- data_ex3 |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  filter(alg%nin%paste0('IIT-RF-',1:5)) |>
  filter(p_1==50) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 3')
ex3 #Show
jpeg(file.path(getwd(),'plots',"ex3.jpg"), width = 700)
ex3
dev.off() 
ex3_2 <- data_ex3 |> 
  mutate(post_calls=ifelse(calls_for_pi>500000,500000,calls_for_pi)) |> 
  filter(alg %in% c(paste0('IIT-RF-',1:5),'IIT')) |>
  filter(p_1==50) |> 
  ggplot(aes(x=theta,y=post_calls)) +
  geom_boxplot(aes(fill=alg))+
  labs(x=TeX("$\\theta$"), y=TeX('Calls for $\\pi$'),title='Example 3')
ex3_2 #Show
jpeg(file.path(getwd(),'plots',"ex3_2.jpg"), width = 700)
ex3_2
dev.off() 


