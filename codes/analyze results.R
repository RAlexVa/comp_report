library(tidyverse)
library(stringr)
##### Consolidating results in a single file #####
data_full <- tibble()
trajectory <- list()
model <- list()
temp_ladder <- list()
for(method in 1:4){
  if(method==1){
    archivos <- files_results[grepl(paste0('^resultados_VT-IIT_modelo',method-1), files_results)]
    temp_consolidado <- tibble()
    for(file in archivos){
      
      start_sim <- gsub("^.*sim(\\d+)_.*$", "\\1", file)
      end_sim <- gsub("(.*_){2}(\\d+)\\.csv$", "\\2", file)
      # print(paste0('Method:',method-1,', sim:',start_sim,'-',end_sim))
      temp_read <- read.csv(paste0('results/',file), header=F) #read file
      #Add labels to identify the 50k iterations of each simulation
      temp <- tibble(sim=rep(start_sim:end_sim, each = 50000),mode=temp_read$V1)
      temp_consolidado <- rbind(temp_consolidado,temp)
    }
    # trajectory[[method]] <- temp_consolidado;
    # model[[method]] <- 'IIT';
    temp_consolidado$temp <- 1
    temp_consolidado$t_ladder <- 0
    temp_consolidado$method <- 0
    data_full <- rbind(data_full,temp_consolidado)
  }else{
    list_files <- files_results[grepl(paste0('^resultados_VT-IIT_modelo',method-1), files_results)]
    temps <- gsub("^.*temp_(\\d+)_.*$", "\\1", list_files)
    for(ttt in unique(temps)){
      archivos <- list_files[temps==ttt]
      temp_consolidado <- tibble()
      for(file in archivos){
        start_sim <- gsub("^.*sim(\\d+)_.*$", "\\1", file)
        end_sim <- gsub("(.*_){2}(\\d+)\\.csv$", "\\2", file)
        # print(paste0('Method:',method-1,', temp:',ttt,', sim:',start_sim,'-',end_sim))
        temp_read <- read.csv(paste0('results/',file), header=F) #read file
        #Add labels to identify the 50k iterations of each simulation
        temp <- tibble(sim=rep(start_sim:end_sim, each = 50000),mode=temp_read$V1, temp=temp_read$V2, t_ladder=ttt)
        temp_consolidado <- rbind(temp_consolidado,temp)
      } 
      temp_consolidado$method <- method;
      data_full <- rbind(data_full,temp_consolidado)
    }
  }
}
dim(data_full)
head(data_full)

saveRDS(data_full,'results/VT_IIT_sim_results.rds')


##### Analyzing results #####

data <- readRDS('results/VT_IIT_sim_results.rds')

### Summarize data
resumen <- data |> select(-temp) |> 
  group_by(method,t_ladder,sim,mode) |>
  summarise(visited=n()) |> ungroup()
resumen$mode <- paste0('m',resumen$mode)

resumen <- resumen |> pivot_wider(names_from = mode,values_from = visited)

resumen[is.na(resumen)] <- 0
# #Check that the row sums coincide
# total_sim <- rowSums(resumen |> select(-method,-t_ladder,-sim))
# max(total_sim);min(total_sim)

percentages <- resumen |> group_by(method,t_ladder,sim) |> 
  mutate(total=sum(c_across(m0:m3))) |> 
  mutate(across(m0:m3,~./total)) |> ungroup()
# Checking how a 1 is presented in percentages
# percentages|> filter(method==4, t_ladder==2, sim==52)

visited <- percentages |> 
           select(-total,-m0)|>
           mutate(across(m4:m3,~ceiling(.))) |> 
           group_by(method,t_ladder,sim) |> 
           mutate(total=sum(c_across(m4:m3)),
                  group1=sum(c_across(m1:m3)),
                  group2=sum(c_across(m4:m6))) |> ungroup()


final_table <- visited |> 
  select(method,t_ladder,total) |> 
  group_by_all() |> 
  summarise(simulations=n()) |> 
  ungroup() |> 
  pivot_wider(names_from = total, values_from = simulations)

final_table[is.na(final_table)] <- 0

final_table$Method <- ifelse(final_table$method==0,'IIT',paste0('M',final_table$method-1))
final_table$Delta <- ceiling(as.numeric(final_table$t_ladder)/2)
final_table$Tot_temps <- 5+5*(1-as.numeric(final_table$t_ladder)%%2)
final_table$Tot_temps <- ifelse(final_table$Method=='IIT',5,final_table$Tot_temps)
final_table <- final_table |> arrange(Delta,Tot_temps,Method)

final_table <- final_table |> 
  select(Method,Delta,Tot_temps,`0`,`1`,`2`,`3`,`4`,`5`,`6`)

final_table |> filter(Tot_temps==5)
final_table$Tot_temps <- ifelse(final_table$Method=='IIT',10,final_table$Tot_temps)
final_table |> filter(Tot_temps==10)
