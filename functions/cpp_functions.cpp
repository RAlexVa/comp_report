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


//Below function works with loglikelihood
// This is the update for method 1
// [[Rcpp::export]]
List VT_IIT_update_c1(colvec X, vec temp, uvec curr_temp, mat modelX, vec resY, int n, int t, vec logpsi) {
  double n_double = double(n);
  double t_double = double(t);
  // Create the model matrix
  double temperature = conv_to<double>::from(temp.elem(curr_temp)); // Current temperature
  // Rcpp::Rcout << temperature << std::endl;
  int total_neighbors = n+t;
  vec probs(total_neighbors, fill::zeros); //log probabilities

// Compute likelihood of the current state
double logpi_current=0;
uvec current_coord = find(X==1);
if(sum(current_coord)==0){ // If the current state is all zeros
  logpi_current=logL_0(resY);
}else{ // If the current state is not all zeroes
  logpi_current = logLikelihood(modelX,resY,current_coord);
}
// Rcpp::Rcout << logpi_current << std::endl;
 
 
//Compute weight for spatial neighbors
double temporal=0;
 for(int j=0; j<n;j++){
   vec newX = X;
   newX.row(j) = 1-X.row(j);
   //Rcpp::Rcout << newX << std::endl;
   uvec coord = find(newX==1);
   //Rcpp::Rcout << coord << std::endl;
   if(sum(coord)==0){
   temporal=(logL_0(resY)-logpi_current)*temperature;
   // Rcpp::Rcout << "all zeroes "<<temporal << std::endl;
   if(temporal<0){probs(j) =exp(temporal)/n_double;}else{probs(j) =1/n_double;}
   }else{
     temporal=((logLikelihood(modelX,resY,coord)-logpi_current)*temperature);
     // Rcpp::Rcout << sum(coord) <<" coordinates "<< temporal << 1/n_double << "que esta pasando" <<std::endl;
     if(temporal<0){probs(j) = exp(temporal)/n_double;}else{probs(j) =1/n_double;}
   }
 }
 
 // Compute weight for temperature neighbors
 double temp_nei = 0;
 for(int j=n; j<n+t;j++){
   temp_nei = conv_to<double>::from(temp.row(j-n)); //Get the value of the temperature
   // Rcpp::Rcout << temp_nei << std::endl;
   //Usually we deleted the entry for the current temperature
   //But it's easier just to assign a probability 0 to choose it
   if(temp_nei==temperature){probs(j)=0;}else{//
     temporal=(logpi_current*(temp_nei-temperature) + conv_to<double>::from(logpsi.row(j-n))-conv_to<double>::from(logpsi.elem(curr_temp)));
      // Rcpp::Rcout << temporal << std::endl;
      // Rcpp::Rcout << logpi_current*(temp_nei-temperature) << std::endl;
     if(temporal<0){probs(j) = (exp(temporal)/(t_double-1));}else{probs(j)=(1/(t_double-1));}}
 }
 // Rcpp::Rcout << probs << std::endl;
 //Choose the next neighbor
 vec u = Rcpp::runif(n+t);
 vec probs_choose = -log(u)/probs;
 // Rcpp::Rcout << u << std::endl;
 // Rcpp::Rcout << log(u) << std::endl;
 // Rcpp::Rcout << -log(u)/probs << std::endl;
 
 //Find the index of the minimum element. source:https://gallery.rcpp.org/articles/vector-minimum/
 //This corresponds to choosing that neighbor
 int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
 // Rcpp::Rcout << neigh_pos << std::endl;
 
 if(neigh_pos<n){
   X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
 }else{
   temperature = neigh_pos-n; //Update the temperature
 }
 
 List ret;
 ret["temp"]=temperature;
 // ret["lprobs"]=probs;
 ret["X"]=X;
 return ret;
}

// [[Rcpp::export]]
vec update_logpsi(vec logpsi, int curr_temp, double iteration, int J){
  
  double n_0=100;
  double s_0=100;
  double temporal = conv_to<double>::from(logpsi.row(curr_temp)) - s_0/(n_0+iteration);
  logpsi = logpsi + s_0/(J*(n_0 + iteration));
  logpsi.row(curr_temp) = temporal;
  return(logpsi);
}


// [[Rcpp::export]]
vec Simulation_mod1(int n, int p, int numsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything
  //Rcpp::Rcout << "Starts function "<< std::endl;
// Repeats according to the number of simulations
for(int s=0; s<numsim;s++){
  //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
  // Creates model
  // Create the model matrix X and response Y
  mat modelX(n,p);
  vec resY(n);
  vec beta;
  // Define beta
  beta = Rcpp::runif(1);
  beta = (4+(beta*2))*sqrt(log(p)/n);
  for(int j=0; j<p;j++){
    modelX.col(j) = as<vec>(Rcpp::rnorm(n));
  }
  // Define response
  resY= ((modelX.col(0) + modelX.col(1) + modelX.col(2))*beta) + as<vec>(Rcpp::rnorm(n));
  
  //Create multimodality
  modelX.col(3) = modelX.col(1)-modelX.col(2)+as<vec>(Rcpp::rnorm(n));
  modelX.col(4) = modelX.col(0)+modelX.col(1)+modelX.col(2)+modelX.col(5)+modelX.col(6)+as<vec>(Rcpp::rnorm(n));
  modelX.col(7) = modelX.col(5)-modelX.col(6)+as<vec>(Rcpp::rnorm(n));
 // Rcpp::Rcout << "Model defined  "<< modelX<<std::endl; 
 // Rcpp::Rcout << "Model defined  "<< resY<<std::endl; 
// Here we have the model defined for iteration #s

///////////////////////////////////////////////////////////////////////////
//Then we start the for loop to run over iterations
vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
// Rcpp::Rcout << "Creates vector of 0  "<< X<<std::endl;
uvec curr_temp = "0"; // We start in the very first temperature (which should be 1)
// Rcpp::Rcout << "Defines initial temperature "<< curr_temp<<std::endl;
vec logpsi(t, fill::zeros); // Initialize vector logpsi, as many entries as temperatures
//Starts simulation
for(int i=0;i<numiter;i++){
  // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 
  double n_double = double(p);
  double t_double = double(t);
  
  // Create the model matrix
  double temperature = conv_to<double>::from(temp.elem(curr_temp)); // Current temperature
  
  int total_neighbors = n+t;
  vec probs(total_neighbors, fill::zeros); //probabilities
  
  // Compute likelihood of the current state
  double logpi_current=0;
  uvec current_coord = find(X==1);
  if(sum(current_coord)==0){ // If the current state is all zeros
    logpi_current=logL_0(resY);
  }else{ // If the current state is not all zeroes
    logpi_current = logLikelihood(modelX,resY,current_coord);
  }
  
  //Compute weight for spatial neighbors
  double temporal=0;
  for(int j=0; j<n;j++){
    Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
    vec newX = X;
    newX.row(j) = 1-X.row(j);
    //Rcpp::Rcout << newX << std::endl;
    uvec coord = find(newX==1);
    //Rcpp::Rcout << coord << std::endl;
    if(sum(coord)==0){
      temporal=(logL_0(resY)-logpi_current)*temperature;
      // Rcpp::Rcout << "all zeroes "<<temporal << std::endl;
      if(temporal<0){probs(j) =exp(temporal)/n_double;}else{probs(j) =1/n_double;}
    }else{
      temporal=((logLikelihood(modelX,resY,coord)-logpi_current)*temperature);
      // Rcpp::Rcout << sum(coord) <<" coordinates "<< temporal << 1/n_double << "que esta pasando" <<std::endl;
      if(temporal<0){probs(j) = exp(temporal)/n_double;}else{probs(j) =1/n_double;}
    }
  }
  
  // Compute weight for temperature neighbors
  double temp_nei = 0;
  for(int j=n; j<n+t;j++){
    temp_nei = conv_to<double>::from(temp.row(j-n)); //Get the value of the temperature
    //Usually we deleted the entry for the current temperature
    //But it's easier just to assign a probability 0 to choose it
    if(temp_nei==temperature){probs(j)=0;}else{//
      temporal=(logpi_current*(temp_nei-temperature));
      if(temporal<0){probs(j) = (exp(temporal)/(t_double-1));}else{probs(j)=(1/(t_double-1));}}
  }
  
  //Choose the next neighbor
  vec u = Rcpp::runif(n+t);
  vec probs_choose = -log(u)/probs;
  
  //Find the index of the minimum element. source:https://gallery.rcpp.org/articles/vector-minimum/
  //This corresponds to choosing that neighbor
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  // Rcpp::Rcout << neigh_pos << std::endl;
  
  if(neigh_pos<n){
    X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  }else{
    temperature = neigh_pos-n; //Update the temperature
  }
  Rcpp::Rcout << "Finish iteration "<< i << " of simulation "<< s << std::endl;
  Rcpp::Rcout << "X= "<< X << std::endl;
  Rcpp::Rcout << "temp= "<< temperature << std::endl;
  
}// End of for loop for iterations
}//End of for loop for simulations

return temp;
}//End of function

// How to use functions inside other functions
//If they're in the same file, they can interact
//https://stackoverflow.com/questions/45296381/can-we-use-rcpp-with-multiple-c-functions/45296773#45296773


//Translating code from R to armadillo



// COnvert vectors from Rcpp to armadillo
// https://statisticsglobe.com/convert-vector-from-rcpparmadillo-to-rcpp-r

