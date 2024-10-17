#Create 100 random models following IIT paper the setup for testing VT-IIT algorithm
library(Rcpp)
library(RcppArmadillo)

## Function to create random models
Rcpp::sourceCpp("functions/cpp_func_random_model.cpp")

set.seed(1110)
numsim <- 100
betas <- numeric(numsim)
for(i in 1:numsim){
  mod <- random_model(n=100,p=200)
  write.table(mod$X,file=paste0('models/modelX',i,'.csv'),sep=',',row.names = F,col.names=F);
  write.table(mod$Y,paste0('models/resY',i,'.csv'),sep=',',row.names = F,col.names=F);
  betas[i] <- mod$beta
}
write.table(betas,file=paste0('models/betas.csv'),sep=',',row.names = F,col.names=F);
#We create this set of matrices and response vectors to use the exact same model when
#testing different algorithms perfmance


#Checking that seed works
# X1v1 <- read.csv('models/modelX1.csv')
# X100v1 <- read.csv('models/modelX100.csv')
# Y1v1 <- read.csv('models/resY1.csv')
# Y100v1 <- read.csv('models/resY100.csv')
# 
# new_path <- 'C:/Users/ralex/Documents/'
# X1v2 <- read.csv(paste0(new_path,'modelX1.csv'))
# X100v2 <- read.csv(paste0(new_path,'modelX100.csv'))
# Y1v2 <- read.csv(paste0(new_path,'resY1.csv'))
# Y100v2 <- read.csv(paste0(new_path,'resY100.csv'))
# identical(X1v1,X1v2)
# identical(X100v1,X100v2)
# identical(Y1v1,Y1v2)
# identical(Y100v1,Y100v2)
