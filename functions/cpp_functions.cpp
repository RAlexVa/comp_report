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
  
//Manual implementation of linear regression
//source: https://cran.r-project.org/web/packages/RcppArmadillo/readme/README.html
int n = X.n_rows;//, k = X.n_cols

arma::colvec coef = arma::solve(subX, Y);     // fit model y ~ X
arma::colvec res  = Y - subX*coef;            // residuals
double rss = arma::dot(res, res);  // 

double logl = -0.5*(n*log(2*M_PI) + n*log(rss) -n*log(n) + n);
return logl;
}
// [[Rcpp::export]]
double logL_0(colvec Y){
  int n = Y.n_rows;
  vec subX(n,fill::ones);
  
  //Manual implementation of linear regression with no features
  arma::colvec coef = arma::solve(subX, Y);     // fit model y ~ X
  arma::colvec res  = Y - subX*coef;            // residuals
  double rss = arma::dot(res, res);  // 
  
  double logl = -0.5*(n*log(2*M_PI) + n*log(rss) -n*log(n) + n);
  return logl;
}


// [[Rcpp::export]]
List VT_IIT_update_c(vec X, vec temp, uvec curr_temp, mat modelM, vec resY,int n, int t) {
  // Create the model matrix
  vec temperature = temp.elem(curr_temp); // Current temperature
  int total_neighbors = n+t;
  vec logprobs(total_neighbors, fill::zeros); //log probabilities
 
//Update vector
 for(int j=0; j<total_neighbors;j++){
   logprobs(j) = j;
 }
 List ret;
 ret["temp"]=temperature;
 ret["lprobs"]=logprobs;
 return ret;

}




/*** R
# VT_IIT_update_c(X=c(1,2),
#                  temp=c(0.5,0.4,0.1),
#                  curr_temp=2,
#                  modelM=matrix(1:4,2,2),
#                  resY=c(1,2),
#                  n=30,
#                  t=5)
*/
