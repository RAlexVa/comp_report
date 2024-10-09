//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

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
//Return the loglikelihood when there are no covariates selected
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

//Update the vector of logpsi
// [[Rcpp::export]]
vec update_logpsi(vec logpsi, int curr_temp, double iteration, int J){
  
  double n_0=100;
  double s_0=100;
  double temporal = conv_to<double>::from(logpsi.row(curr_temp)) - s_0/(n_0+iteration);
  logpsi = logpsi + s_0/(J*(n_0 + iteration));
  logpsi.row(curr_temp) = temporal;
  return(logpsi);
}

// Function to identify if one of the 6 modes were visited
// [[Rcpp::export]]
int check_modes(const arma::vec& X) {
  // Defined modes
  arma::vec mod1 = {0, 1, 2};
  arma::vec mod2 = {0, 1, 3};
  arma::vec mod3 = {0, 2, 3};
  arma::vec mod4 = {4, 5, 6};
  arma::vec mod5 = {4, 5, 7};
  arma::vec mod6 = {4, 6, 7};
  
  // Collect vectors into a list
  std::vector<arma::vec> modes = {mod1, mod2, mod3, mod4, mod5, mod6};
  
  vec pos=conv_to<vec>::from(find(X==1));
  for (size_t i = 0; i < modes.size(); i++) {
    //Rcpp::Rcout << modes[i]<< std::endl;
    if (arma::approx_equal(pos, modes[i], "absdiff", 0)) {  // Exact comparison
      return i; // Return true if X is equal to any model
    }
  }
  return -1; // Return false if X does not match any model
}

// [[Rcpp::export]]
// This function keeps returns all the trajectory for the process considering all iterations
//But only considers the latest simulation, for every new s the output is refreshed
List Simulation_mod1_full(int n, int p, int numsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything
  //Rcpp::Rcout << "Starts function "<< std::endl;
//////////////////////////////////////////////////
// Define modes to check how many times they're visited
//modes are: 1,2,3;  1,2,4; 1,3,4; 5,6,7; 5,6,8; 5,7,8;


  //To store the results
  mat states_visited(numiter,p);//Matrix to store the states visited
  vec temps_visited(numiter);//Vector to store the temperatures visited
  mat logpsi_record(numiter,t);//Vector to store the values of logpsi used  
// Repeats according to the number of simulations
/////////////////////////////////////////////////////
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
//Reset the variables to store the results
states_visited.zeros();
temps_visited.zeros();
logpsi_record.zeros();
//Then we start the for loop to run over iterations
vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
// Rcpp::Rcout << "Creates vector of 0  "<< X<<std::endl;
int curr_temp = 0; // We start in the very first temperature (which should be 1)
// Rcpp::Rcout << "Defines initial temperature "<< curr_temp<<std::endl;
vec logpsi(t, fill::zeros); // Initialize vector logpsi, as many entries as temperatures
//Starts simulation
for(int i=0;i<numiter;i++){
  // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 
  double n_double = double(p);
  double t_double = double(t);
  
  // Create the model matrix
  double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
  
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
    // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
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
      // temporal=(logpi_current*(temp_nei-temperature)); //This was without the logpsi factors
      temporal=(logpi_current*(temp_nei-temperature) + conv_to<double>::from(logpsi.row(j-n))-conv_to<double>::from(logpsi.row(curr_temp)));
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
    curr_temp = neigh_pos-n; //Update the temperature
  }
  
// Update logpsi vector
  logpsi = update_logpsi(logpsi,curr_temp,double(i),t-1);
  // Rcpp::Rcout << logpsi << std::endl;
  
  //update_logpsi(vec logpsi, int curr_temp, double iteration, int J)
  // Rcpp::Rcout << "Finish iteration "<< i << " of simulation "<< s << std::endl;
  // Rcpp::Rcout << "X= "<< X << std::endl;
  // Rcpp::Rcout << "temp= "<< temperature << std::endl;
  
  
// Store the new values to report
states_visited.row(i)=X.t();//Matrix to store the states visited
temps_visited(i)=curr_temp;//Vector to store the temperatures visited
logpsi_record.row(i)=logpsi.t();//Vector to store the values of logpsi used
}// End of for loop for iterations
// Rcpp::Rcout << "Finish simulation "<< s << std::endl;
// Rcpp::Rcout << "X= "<< X << std::endl;
}//End of for loop for simulations

List ret;
  ret["states"]=states_visited;
  ret["temps"]=temps_visited;
  ret["logpsi"]=logpsi_record;
  return ret;

}//End of function

// [[Rcpp::export]]
// This function keeps returns all the trajectory for the process considering all iterations
//But only considers the latest simulation, for every new s the output is refreshed
vec Simulation_mod1(int n, int p, int numsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything
  //Rcpp::Rcout << "Starts function "<< std::endl;
  //////////////////////////////////////////////////
  // Define modes to check how many times they're visited
  //modes are: 1,2,3;  1,2,4; 1,3,4; 5,6,7; 5,6,8; 5,7,8;
  
  
  //To store the results
  vec modes_visited(numiter * numsim);//Vector to store the modes visited
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
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
    int curr_temp = 0; // We start in the very first temperature (which should be 1)
    // Rcpp::Rcout << "Defines initial temperature "<< curr_temp<<std::endl;
    vec logpsi(t, fill::zeros); // Initialize vector logpsi, as many entries as temperatures
    //Starts simulation
    for(int i=0;i<numiter;i++){
      if (i % 1000 == 1) {  // equivalent to i%%1000==1 in R
        
        Rcpp::Rcout << "Simulation: " << s << "Iteration: " << i << std::endl;
      }
      // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 
      double n_double = double(p);
      double t_double = double(t);
      
      // Create the model matrix
      double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
      
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
        // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
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
          // temporal=(logpi_current*(temp_nei-temperature)); //This was without the logpsi factors
          temporal=(logpi_current*(temp_nei-temperature) + conv_to<double>::from(logpsi.row(j-n))-conv_to<double>::from(logpsi.row(curr_temp)));
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
        curr_temp = neigh_pos-n; //Update the temperature
      }
      
      // Update logpsi vector
      logpsi = update_logpsi(logpsi,curr_temp,double(i),t-1);
      // Rcpp::Rcout << logpsi << std::endl;
      
      //update_logpsi(vec logpsi, int curr_temp, double iteration, int J)
      // Rcpp::Rcout << "Finish iteration "<< i << " of simulation "<< s << std::endl;
      // Rcpp::Rcout << "X= "<< X << std::endl;
      // Rcpp::Rcout << "temp= "<< temperature << std::endl;
      
      
      // Store the new values to report
      modes_visited.row((s*numiter)+i)=check_modes(X)+1;
    }// End of for loop for iterations
    // Rcpp::Rcout << "Finish simulation "<< s << std::endl;
    // Rcpp::Rcout << "X= "<< X << std::endl;
  }//End of for loop for simulations

  return modes_visited;
}//End of function

