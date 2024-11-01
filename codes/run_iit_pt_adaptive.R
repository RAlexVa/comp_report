#if(!require('Rcpp')){install.packages('Rcpp')}
library(Rcpp)
library(RcppArmadillo)
setwd('..')
Rcpp::sourceCpp("functions/cpp_functions.cpp")

#Define parameters for the model
seed_def <- 6055
# n <- 100
p <- 200
# iterations <- 20000

writeLines('1 es Adaptive\n2 es bounded')
algorithm <- as.numeric(readline('Select algorithm'))

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

writeLines('method 1 min bal. fun. for all\nmethod 2 half & half bal. fun.\nmethod 3 only max temp is different')
m_selected <- paste0("M",as.numeric(readline('Select method')))

for(chunk_selected in 1:5){
  start_point <- 1+(chunk_selected-1)*20
  end_point <- chunk_selected*20
  set.seed(seed_def + chunk_selected)
  if(algorithm==1){
    results <- PT_IIT_adapt_sim(p= p,startsim=start_point,endsim=end_point,L_samples=1000,total_swaps=300, temp=temp, m_selected)
    write.table(results$modes,paste0('results/','resultados_ada_PT-IIT_modelo',m_selected,'_temp_',t_selected,'_seed_',seed_def,'+',chunk_selected,'sim',start_point,'_',end_point,'_modes','.csv'),row.names=F, col.names=F, sep=',')  
    write.table(results$ip,paste0('results/','resultados_ada_PT-IIT_modelo',m_selected,'_temp_',t_selected,'_seed_',seed_def,'+',chunk_selected,'sim',start_point,'_',end_point,'_ip','.csv'),row.names=F, col.names=F, sep=',')   
  }
  if(algorithm==2){
    results <- PT_IIT_bounded_sim(p= p,startsim=start_point,endsim=end_point,L_samples=1000,total_swaps=300, temp=temp, m_selected, prob_logbound = rep(4,length(temp)))    
    write.table(results$modes,paste0('results/','resultados_bound_PT-IIT_modelo',m_selected,'_temp_',t_selected,'_seed_',seed_def,'+',chunk_selected,'sim',start_point,'_',end_point,'_modes','.csv'),row.names=F, col.names=F, sep=',')  
    write.table(results$ip,paste0('results/','resultados_bound_PT-IIT_modelo',m_selected,'_temp_',t_selected,'_seed_',seed_def,'+',chunk_selected,'sim',start_point,'_',end_point,'_ip','.csv'),row.names=F, col.names=F, sep=',')   
  }

}



