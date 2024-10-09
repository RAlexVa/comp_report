if(!require('Rcpp')){install.packages('Rcpp')}
library(Rcpp)
setwd('..')
Rcpp::sourceCpp("functions/cpp_functions.cpp")

n <- 100
p <- 200
iterations <- 50000
simulations <- 50
temp1.1 <- (1+((1:5)-1))
temp <- temp1.1
results <- Simulation_mod1(n=n,p=p,numsim=simulations,numiter=iterations,temp=temp,t=length(temp))

write.csv(results,'results/resultados_modelo1.csv',row.names=F)
