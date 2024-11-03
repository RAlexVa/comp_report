library(tidyverse)
library(stringr)
##### For VT-IIT #####
##### Consolidating results in a single file #####
data_full <- tibble()
files_results <- dir('results')
# check <- files_results[grepl('^resultados_VT-IIT_modelo\\d', files_results)]
for(method in 0:3){
  print(paste0("method: ",method))
  if(method==0){
    archivos <- files_results[grepl(paste0('^resultados_VT-IIT_modelo',method), files_results)]
    temp_consolidado <- tibble()
    for(file in archivos){
      start_sim <- as.numeric(gsub("^.*sim(\\d+)_.*$", "\\1", file))
      end_sim <- as.numeric(gsub("^.*sim(\\d+)_(\\d+).*$", "\\2", file))
      if(grepl("it_", file)){
        numit <- as.numeric(gsub(".*it_([0-9]+)\\.csv.*", "\\1", file))
      }else{numit <- -1}
      # print(paste0('Method:',method-1,', sim:',start_sim,'-',end_sim))
      temp_read <- read.csv(paste0('results/',file), header=F) #read file
      #Add labels to identify the 50k iterations of each simulation
      temp <- tibble(sim=rep(start_sim:end_sim, each = 50000*abs(numit)),mode=temp_read$V1)
      temp$temp <- 1
      temp$t_ladder <- 0
      temp$method <- 0
      temp$iterations <- numit
      temp_consolidado <- rbind(temp_consolidado,temp)
    }
    saveRDS(temp_consolidado,'results/VT_IIT_sim_results_m0.rds')
    # data_full <- rbind(data_full,temp_consolidado)
  }else{
    list_files <- files_results[grepl(paste0('^resultados_VT-IIT_modelo',method), files_results)]
    temps <- gsub("^.*temp_(\\d+)_.*$", "\\1", list_files)
    temp_consolidado <- tibble()
    for(ttt in unique(temps)){
      print(paste0("temp: ",ttt))
      archivos <- list_files[temps==ttt]
      for(file in archivos){
        print(paste0("file: ",file))
        start_sim <- as.numeric(gsub("^.*sim(\\d+)_.*$", "\\1", file))
        end_sim <- as.numeric(gsub("^.*sim(\\d+)_(\\d+).*$", "\\2", file))
        if(grepl("it_", file)){
          numit <- as.numeric(gsub(".*it_([0-9]+)\\.csv.*", "\\1", file))
        }else{numit <- -1}
        # print(paste0('Method:',method-1,', temp:',ttt,', sim:',start_sim,'-',end_sim))
        temp_read <- read.csv(paste0('results/',file), header=F) #read file
        #Add labels to identify the 50k iterations of each simulation
        temp <- tibble(sim=rep(start_sim:end_sim, each = 50000*abs(numit)),mode=temp_read$V1, temp=temp_read$V2, t_ladder=ttt)
        temp$method <- method;
        temp$iterations <- numit
        temp_consolidado <- rbind(temp_consolidado,temp)
      } 
      # data_full <- rbind(data_full,temp_consolidado)
    }
    saveRDS(temp_consolidado,paste0('results/VT_IIT_sim_results_m',method,'.rds'))
  }
}
# dim(data_full)
# head(data_full)
# 
# saveRDS(data_full,'results/VT_IIT_sim_results.rds')


##### Analyzing results #####

# data <- readRDS('results/VT_IIT_sim_results.rds')
res_tot <- NA
for(method in 0:3){
  data <- readRDS(paste0('results/VT_IIT_sim_results_m',method,'.rds'))
  ### Summarize data
  resumen <- data |> select(-temp) |> 
    group_by(method,t_ladder,sim,mode, iterations) |>
    summarise(visited=n()) |> ungroup()
  resumen$mode <- paste0('m',resumen$mode)
  
  resumen <- resumen |> pivot_wider(names_from = mode,values_from = visited)
  res_tot <- rbind(res_tot,resumen)
}
resumen <- res_tot[-1,]

rm(list=c('res_tot'))
resumen[is.na(resumen)] <- 0
#dim(resumen)
#100*3*4*2 + 100*2
# #Check that the row sums coincide
 # total_sim <- rowSums(resumen |> select(-method,-t_ladder,-sim,-iterations))
 # max(total_sim);min(total_sim)

percentages <- resumen |> group_by(method,t_ladder,sim,iterations) |> 
  mutate(total=sum(c_across(m0:m3))) |> 
  mutate(across(m0:m3,~./total)) |> ungroup()
# Checking how a 1 is presented in percentages
# percentages|> filter(method==4, t_ladder==2, sim==52)

visited <- percentages |> 
           select(-total,-m0)|>
           mutate(across(m4:m3,~ceiling(.))) |> 
           group_by(method,t_ladder,sim,iterations) |> 
           mutate(total=sum(c_across(m4:m3)),
                  group1=sum(c_across(m1:m3)),
                  group2=sum(c_across(m4:m6))) |> ungroup()

final_table <- visited |> 
  select(method,t_ladder,iterations,total) |> 
  group_by_all() |> 
  summarise(simulations=n()) |> 
  ungroup() |> 
  pivot_wider(names_from = total, values_from = simulations)

# final_table <- visited |> select(-group1,-group2,-total) |> 
#   group_by(method,t_ladder,iterations) |> 
#   summarise(across(starts_with("m"), \(x)sum(x, na.rm=T))) |> 
#   ungroup() 

final_table[is.na(final_table)] <- 0

final_table$Method <- ifelse(final_table$method==0,'IIT',paste0('M',final_table$method))
final_table$Delta <- ceiling(as.numeric(final_table$t_ladder)/2)
final_table$Tot_temps <- 5+5*(1-as.numeric(final_table$t_ladder)%%2)
final_table$Tot_temps <- ifelse(final_table$Method=='IIT',5,final_table$Tot_temps)
final_table <- final_table |> arrange(Delta,Tot_temps,Method)

final_table <- final_table |> arrange(iterations) |> 
  select(Method,Delta,Tot_temps,iterations,`1`,`2`,`3`,`4`,`5`,`6`)

final_table |> filter(Tot_temps==5)
final_table |> filter(Tot_temps==10|Method=='IIT')


##### Exporting for VT-IIT #####
saveRDS(final_table,paste0('results/VT_IIT_full_results.rds'))


vtiit <- readRDS(paste0('results/VT_IIT_full_results.rds'))




##### For PT-IIT #####
library(stringr)
data_full <- tibble()
files_results <- dir('results')
files_results <- files_results[grepl("PT-IIT",files_results)]
files_results <- files_results[grepl("_modes\\.csv$",files_results)]

for(i in 1:length(files_results)){
  file <- files_results[i]#1,40,80,120
  alg <- str_extract(file, "(?<=resultados_).*?(?=_PT-IIT)")
  start_sim <- as.numeric(gsub("^.*sim(\\d+)_.*$", "\\1", file))
  end_sim <- as.numeric(gsub("^.*sim(\\d+)_(\\d+).*$", "\\2", file))
  model <- str_extract(file, "(?<=modelo).{2}")
  temperature <- str_extract(file, "(?<=temp_).{1}")
  temp_read <- read.csv(paste0('results/',file), header=F)
  
  
  
  if(is.na(alg)|alg=='noad'){
    if(is.na(alg)){alg <- 'pt'}
    if(nrow(temp_read)!=400000){print(paste0("Wrong dimensions in ",i))}
    temp_read$sim <- rep(start_sim:end_sim, each = 20000)
    
    temp_read <- temp_read |> 
      pivot_longer(-sim, names_to='replica',values_to = 'mode') |> 
      select(-replica) |> 
      group_by(sim,mode) |> 
      summarise(visits=n()) |> 
      ungroup() |> 
      pivot_wider(id_cols=sim,names_from=mode, values_from=visits)
    
    if(ncol(temp_read)!=8){print(paste0("Less columns in ",i))}
    temp_read[is.na(check)] <- 0
    colnames(temp_read) <- c("sim",paste0('m',0:6))
    temp_read$model <- model
    temp_read$alg <- alg
    temp_read$temp <- temperature
  }
  
  if(alg %in% c('ada','bound')){
    if(nrow(temp_read)!=6000){print(paste0("Wrong dimensions in ",i))}
    temp_read$sim <- rep(start_sim:end_sim, each = 300)
    
    temp_read <- temp_read |> group_by(sim) |>
      summarise(across(everything(), \(x) sum(x, na.rm = TRUE)))
    if(ncol(temp_read)!=8){print(paste0("Less columns in ",i))}
    colnames(temp_read) <- c("sim",paste0('m',0:6))
    temp_read$model <- model
    temp_read$alg <- alg
    temp_read$temp <- temperature
  }
  
  data_full <- rbind(data_full,temp_read)
  
}
data_full[is.na(data_full)] <- 0
saveRDS(data_full,paste0('results/PT_results.rds'))
write.csv(data_full,paste0('results/PT_results.csv'), row.names = F)
##### Checking speed of finding the modes #####
# I need to find the smallest row index
#But row indexes only go from 1 to 50k or 100k depending on iterations
res_tot <- NA
for(method in 0:3){
  data <- readRDS(paste0('results/VT_IIT_sim_results_m',method,'.rds'))
  ### Summarize data
  resumen <- data |> select(-temp) |>  
    mutate(r_i=row_number()) |> 
    # group_by(method,t_ladder,mode,sim,iterations) |> 
    group_by(method,t_ladder,sim,iterations) |> 
    slice_min(r_i,n=1) |> ungroup()
}