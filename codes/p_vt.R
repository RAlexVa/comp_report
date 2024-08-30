rm(list=ls())
setwd('..')

m <- as.numeric(readline('Select method'))


if(m==1){
  source(file.path(getwd(),'codes','example_multimodal.R'))
  ## Method 1
  method <- 1
  for(J in c(4,5,9,10)){
    #Create vector of temperatures depending on J
    logh_func <- list()
    for (j in 1:J)
    {logh_func[[j]] <- hsq_log}
    #Run simulations for the deltas
    for(delta in 1:2){
      set.seed(453)
      VT_IIT(method=method,h_func=logh_func,J=J,delta=delta) 
    }
  }  
}

if(m==2){
  ##########################################################
  source(file.path(getwd(),'codes','example_multimodal.R'))
  ## Method 2
  method <- 2
  for(J in c(4,5,9,10)){
    #Create vector of temperatures depending on J
    logh_func <- list()
    for (j in 1:floor(J/2))
    {logh_func[[j]] <- hsq_log}
    for (j in (floor(J/2)+1):J)
    {logh_func[[j]] <- hmin_log}
    #Run simulations for the deltas
    for(delta in 1:2){
      set.seed(453)
      VT_IIT(method=method,h_func=logh_func,J=J,delta=delta) 
    }}  
}

if(m==3){
  ##########################################################
  source(file.path(getwd(),'codes','example_multimodal.R'))
  ## Method 3
  method <- 3
  for(J in c(4,5,9,10)){
    #Create vector of temperatures depending on J
    logh_func <- list()
    for (j in 1:J)
    {logh_func[[j]] <- hsq_log}
    logh_func[[J]] <- hmin_log
    #Run simulations for the deltas
    for(delta in 1:2){
      set.seed(453)
      VT_IIT(method=method,h_func=logh_func,J=J,delta=delta) 
    }} 
}

if(m==4){
  ##########################################################
  source(file.path(getwd(),'codes','example_multimodal.R'))
  #Method IIT
  method <- 'IIT'
  for(J in c(4,5,9,10)){
    #Create vector of temperatures depending on J
    #Run simulations for the deltas
    for(delta in 1:2){
      set.seed(453)
      just_IIT(method=method,h_func=hsq_log,J=J,delta=delta)
    }} 
}

