library(tidyverse)
ar <- read.csv('results/resultados_VT-IIT_modelo0_seed_6055+1sim1_20.csv', header = F)
data <- tibble(model=sort(rep(1:20,50000)),modes=ar[,1])
