//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//Given dimensions, return a random model 
//following the description provided in the IIT paper
// [[Rcpp::export]]
List random_model(int n, int p) {
// Create the model matrix
  mat ranmat(n,p);
  vec Y(n);
  vec beta;
  // Define beta
  beta = Rcpp::runif(1);
  beta = (4+(beta*2))*sqrt(log(p)/n);
  for(int j=0; j<p;j++){
    ranmat.col(j) = as<vec>(Rcpp::rnorm(n));
  }
  // Define response
  Y= ((ranmat.col(0) + ranmat.col(1) + ranmat.col(2))*beta) + as<vec>(Rcpp::rnorm(n));
  
  //Create multimodality
  
  ranmat.col(3) = ranmat.col(1)-ranmat.col(2)+as<vec>(Rcpp::rnorm(n));
  ranmat.col(4) = ranmat.col(0)+ranmat.col(1)+ranmat.col(2)+ranmat.col(5)+ranmat.col(6)+as<vec>(Rcpp::rnorm(n));
  ranmat.col(7) = ranmat.col(5)-ranmat.col(6)+as<vec>(Rcpp::rnorm(n));
  
  List ret;
  ret["Y"]=Y;
  ret["X"]=ranmat;
  return ret;
}


// Given a model matrix, a response vector and a vector of positions
//return the likelihood of the model using the specified variables
// [[Rcpp::export]]
double logLikelihood(mat X, colvec Y, uvec pos){
  mat subX = X.cols(pos);
  
//Manual implementation of linear regresion
//source: https://cran.r-project.org/web/packages/RcppArmadillo/readme/README.html
int n = X.n_rows;//, k = X.n_cols

arma::colvec coef = arma::solve(subX, Y);     // fit model y ~ X
arma::colvec res  = Y - subX*coef;            // residuals
double rss = arma::dot(res, res);  // 

double logl = -0.5*(n*log(2*M_PI) + n*log(rss) -n*log(n) + n);
return logl;
}

/*** R
#logLikelihood(res$X,res$Y,c(1:3,9:19)-1)
# set.seed(123)
# m <- random_norm(5,4)
# set.seed(123)
# tod <- random_model(4,8)
# 
# r_mod_r <- function(n,p){
#   ##### Defining model #####
#   L <- matrix(rnorm(n*p),nrow=n,ncol=p)
#   
#   beta <- runif(1,4,6)*sqrt(log(p)/n)
#   Y <- (L[,1]+L[,2]+L[,3])*beta + rnorm(n,mean=0,sd=0.5)
#   # Introducing multi modality
#   L[,4] <- L[,2]-L[,3]+rnorm(n,mean=0,sd=0.1)
#   L[,5] <- L[,1]+L[,2]+L[,3]+L[,6]+L[,7]+rnorm(n,mean=0,sd=0.1)
#   L[,8] <- L[,6]-L[,7]+rnorm(n,mean=0,sd=0.1)
#   return(L)
# }
# 
# microbenchmark(random_full(100,200),r_mod_r(100,200),times=200)

*/
