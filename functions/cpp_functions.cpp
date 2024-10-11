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
// This function reports which modes were visited in which iteration and the temperature
// This is for model 1 which uses min balancing function for all temperatures
// [[Rcpp::export]]
mat Simulation_mod1(int n, int p, int numsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything

  //////////////////////////////////////////////////
  mat modes_visited(numiter * numsim,2);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  double t_double = double(t);//
  // Rcpp::Rcout << "t= "<< t <<std::endl;
  // Rcpp::Rcout << "t_double= "<< t_double <<std::endl;
  double J=t_double-1;//Number of temperatures minus 1
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  for(int s=0; s<numsim;s++){
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
//// Creates model
    // Create the model matrix X and response Y
    mat modelX(n,p); // n rows and p columns (p corresponds to the dimension of the problem)
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
    
    //Create multi modality
    modelX.col(3) = modelX.col(1)-modelX.col(2)+as<vec>(Rcpp::rnorm(n));
    modelX.col(4) = modelX.col(0)+modelX.col(1)+modelX.col(2)+modelX.col(5)+modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    modelX.col(7) = modelX.col(5)-modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    // Rcpp::Rcout << "Model defined  "<< modelX<<std::endl; 
    // Rcpp::Rcout << "Model defined  "<< resY<<std::endl; 
    // Here we have the model defined for iteration #s
    
    ///////////////////////////////////////////////////////////////////////////
    //Then we start the for loop to run over iterations
    vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
    // X.row(0)=1;
    // X.row(1)=1;
    // X.row(2)=1;
     // Rcpp::Rcout << "Creates initial vector  "<< X<<std::endl;
     // Rcpp::Rcout << "Check modes initial  "<< check_modes(X)<<std::endl;
    int curr_temp = 0; // We start in the very first temperature (which should be 1)
    // Rcpp::Rcout << "Defines initial temperature "<< curr_temp<<std::endl;
    vec logpsi(t, fill::zeros); // Initialize vector logpsi, as many entries as temperatures
    //Starts simulation
    for(int i=0;i<numiter;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 


      double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
      
      int total_neighbors = p+t; // total number of neighbors is p spacial + t temperature
      // Rcpp::Rcout << "p= "<< p <<std::endl;
      // Rcpp::Rcout << "t= "<< t <<std::endl;
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(sum(current_coord)==0){ // If the current state is all zeros
        logpi_current=logL_0(resY);
      }else{ // If the current state is not all zeroes
        logpi_current = logLikelihood(modelX,resY,current_coord);
      }
////////////      
      //Compute weight for spatial neighbors
      double temporal=0;
      vec newX;
      for(int j=0; j<p;j++){
         // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
        newX = X;
        newX.row(j) = 1-X.row(j);
        //Rcpp::Rcout << newX << std::endl;
        uvec coord = find(newX==1);
        //Rcpp::Rcout << coord << std::endl;
        if(sum(coord)==0){// In case the state visited is all 0s
          temporal=(logL_0(resY)-logpi_current)*temperature;
        }else{// For every other state that is not all 0s
          temporal=((logLikelihood(modelX,resY,coord)-logpi_current)*temperature);
        }
//Apply balancing function to spatial neighbors /////////////////////////////        
if(temporal<0){probs(j) = temporal-log(p_double);}else{probs(j) = -log(p_double);} 
      }
////////////      
      // Compute weight for temperature neighbors
      double temp_nei = 0;
      for(int j=p; j<(p+t);j++){
        // Rcpp::Rcout << "Checks temp neighbors  "<< j<<std::endl; 
        // Rcpp::Rcout << "Checks next neighbor  "<< probs(j+1)<<std::endl; 
        // Rcpp::Rcout << "index for temp  "<< j-p<<std::endl; 
        temp_nei = conv_to<double>::from(temp.row(j-p)); //Get the value of the temperature
        // Rcpp::Rcout << "temperature  "<< temp.row(j-p)<<std::endl;
        //Usually we deleted the entry for the current temperature
        //But it's easier just to assign a probability 0 to choose it
        if(temp_nei==temperature){probs(j)=0;}else{//
          // temporal=(logpi_current*(temp_nei-temperature)); //This was without the logpsi factors
          temporal=(logpi_current*(temp_nei-temperature) + conv_to<double>::from(logpsi.row(j-p))-conv_to<double>::from(logpsi.row(curr_temp)));
//Apply balancing function to temperature neighbors /////////////////////////////           
          if(temporal<0){probs(j) = temporal - log(J);}else{probs(j)=-log(J);}
          }
        // Rcpp::Rcout << "Finish iteration of loop  "<< j<<std::endl; 
      }
      
      //Choose the next neighbor
      vec u = Rcpp::runif(p+t);
      vec probs_choose = log(-log(u))-probs;
      
      //Find the index of the minimum element. source:https://gallery.rcpp.org/articles/vector-minimum/
      //This corresponds to choosing that neighbor
      int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
      // Rcpp::Rcout << neigh_pos << std::endl;
      
      if(neigh_pos<p){
        X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
      }else{
        curr_temp = neigh_pos-p; //Update the temperature
      }
      
      // Update logpsi vector
      logpsi = update_logpsi(logpsi,curr_temp,double(i),t-1);
      // Rcpp::Rcout << logpsi << std::endl;
      //update_logpsi(vec logpsi, int curr_temp, double iteration, int J)

      // Store the new values to report
      // modes_visited.row((s*numiter)+i)=check_modes(X)+1;
      modes_visited.submat((s*numiter)+i,0,(s*numiter)+i,0)=check_modes(X)+1;
      modes_visited.submat((s*numiter)+i,1,(s*numiter)+i,1)=curr_temp+1;
    }// End of for loop for iterations
    // Rcpp::Rcout << "Finish simulation "<< s << std::endl;
    // Rcpp::Rcout << "X= "<< X << std::endl;
  }//End of for loop for simulations
  return modes_visited;
}//End of function

// This function reports which modes were visited in which iteration and the temperature
// This is for model 2 which uses min and sq balancing function half and half of temperatures
// [[Rcpp::export]]
mat Simulation_mod2(int n, int p, int numsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything
  
  //////////////////////////////////////////////////
  mat modes_visited(numiter * numsim,2);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  double t_double = double(t);//
  // Rcpp::Rcout << "t= "<< t <<std::endl;
  // Rcpp::Rcout << "t_double= "<< t_double <<std::endl;
  double J=t_double-1;//Number of temperatures minus 1
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  for(int s=0; s<numsim;s++){
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
    //// Creates model
    // Create the model matrix X and response Y
    mat modelX(n,p); // n rows and p columns (p corresponds to the dimension of the problem)
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
    
    //Create multi modality
    modelX.col(3) = modelX.col(1)-modelX.col(2)+as<vec>(Rcpp::rnorm(n));
    modelX.col(4) = modelX.col(0)+modelX.col(1)+modelX.col(2)+modelX.col(5)+modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    modelX.col(7) = modelX.col(5)-modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    // Rcpp::Rcout << "Model defined  "<< modelX<<std::endl; 
    // Rcpp::Rcout << "Model defined  "<< resY<<std::endl; 
    // Here we have the model defined for iteration #s
    
    ///////////////////////////////////////////////////////////////////////////
    //Then we start the for loop to run over iterations
    vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
    // X.row(0)=1;
    // X.row(1)=1;
    // X.row(2)=1;
    // Rcpp::Rcout << "Creates initial vector  "<< X<<std::endl;
    // Rcpp::Rcout << "Check modes initial  "<< check_modes(X)<<std::endl;
    int curr_temp = 0; // We start in the very first temperature (which should be 1)
    // Rcpp::Rcout << "Defines initial temperature "<< curr_temp<<std::endl;
    vec logpsi(t, fill::zeros); // Initialize vector logpsi, as many entries as temperatures
    //Starts simulation
    for(int i=0;i<numiter;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 
      
      
      double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
      
      int total_neighbors = p+t; // total number of neighbors is p spacial + t temperature
      // Rcpp::Rcout << "p= "<< p <<std::endl;
      // Rcpp::Rcout << "t= "<< t <<std::endl;
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(sum(current_coord)==0){ // If the current state is all zeros
        logpi_current=logL_0(resY);
      }else{ // If the current state is not all zeroes
        logpi_current = logLikelihood(modelX,resY,current_coord);
      }
      ////////////      
      //Compute weight for spatial neighbors
      double temporal=0;
      vec newX;
      for(int j=0; j<p;j++){
        // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
        newX = X;
        newX.row(j) = 1-X.row(j);
        //Rcpp::Rcout << newX << std::endl;
        uvec coord = find(newX==1);
        //Rcpp::Rcout << coord << std::endl;
        if(sum(coord)==0){// In case the state visited is all 0s
          temporal=(logL_0(resY)-logpi_current)*temperature;
        }else{// For every other state that is not all 0s
          temporal=((logLikelihood(modelX,resY,coord)-logpi_current)*temperature);
        }
//Apply balancing function to spatial neighbors /////////////////////////////
if(curr_temp<J/2){//For the first half of temperatures apply the sq balancing function
  probs(j)=(temporal/2) - log(p_double);
}else{//for the other half apply the min balancing function
  if(temporal<0){probs(j) = temporal-log(p_double);}else{probs(j) =-log(p_double);}   
}

      }
      ////////////      
      // Compute weight for temperature neighbors
      double temp_nei = 0;
      for(int j=p; j<(p+t);j++){
        // Rcpp::Rcout << "Checks temp neighbors  "<< j<<std::endl; 
        // Rcpp::Rcout << "Checks next neighbor  "<< probs(j+1)<<std::endl; 
        // Rcpp::Rcout << "index for temp  "<< j-p<<std::endl; 
        temp_nei = conv_to<double>::from(temp.row(j-p)); //Get the value of the temperature
        // Rcpp::Rcout << "temperature  "<< temp.row(j-p)<<std::endl;
        //Usually we deleted the entry for the current temperature
        //But it's easier just to assign a probability 0 to choose it
        if(temp_nei==temperature){probs(j)=0;}else{//
          // temporal=(logpi_current*(temp_nei-temperature)); //This was without the logpsi factors
          temporal=(logpi_current*(temp_nei-temperature) + conv_to<double>::from(logpsi.row(j-p))-conv_to<double>::from(logpsi.row(curr_temp)));
          //Apply balancing function to temperature neighbors /////////////////////////////           
          if(temporal<0){probs(j) = temporal-log(J);}else{probs(j)=-log(J);}
        }
        // Rcpp::Rcout << "Finish iteration of loop  "<< j<<std::endl; 
      }
      
      //Choose the next neighbor
      vec u = Rcpp::runif(p+t);
      vec probs_choose = log(-log(u)) - probs;
      
      //Find the index of the minimum element. source:https://gallery.rcpp.org/articles/vector-minimum/
      //This corresponds to choosing that neighbor
      int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
      // Rcpp::Rcout << neigh_pos << std::endl;
      
      if(neigh_pos<p){
        X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
      }else{
        curr_temp = neigh_pos-p; //Update the temperature
      }
      
      // Update logpsi vector
      logpsi = update_logpsi(logpsi,curr_temp,double(i),t-1);
      // Rcpp::Rcout << logpsi << std::endl;
      //update_logpsi(vec logpsi, int curr_temp, double iteration, int J)
      
      // Store the new values to report
      // modes_visited.row((s*numiter)+i)=check_modes(X)+1;
      modes_visited.submat((s*numiter)+i,0,(s*numiter)+i,0)=check_modes(X)+1;
      modes_visited.submat((s*numiter)+i,1,(s*numiter)+i,1)=curr_temp+1;
    }// End of for loop for iterations
    // Rcpp::Rcout << "Finish simulation "<< s << std::endl;
    // Rcpp::Rcout << "X= "<< X << std::endl;
  }//End of for loop for simulations
  return modes_visited;
}//End of function

// This function reports which modes were visited in which iteration and the temperature
// This is for model 3 which uses min only in temperature J and sq balancing function for every other
// [[Rcpp::export]]
mat Simulation_mod3(int n, int p, int numsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything
  
  //////////////////////////////////////////////////
  mat modes_visited(numiter * numsim,2);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  double t_double = double(t);//
  // Rcpp::Rcout << "t= "<< t <<std::endl;
  // Rcpp::Rcout << "t_double= "<< t_double <<std::endl;
  double J=t_double-1;//Number of temperatures minus 1
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  for(int s=0; s<numsim;s++){
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
    //// Creates model
    // Create the model matrix X and response Y
    mat modelX(n,p); // n rows and p columns (p corresponds to the dimension of the problem)
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
    
    //Create multi modality
    modelX.col(3) = modelX.col(1)-modelX.col(2)+as<vec>(Rcpp::rnorm(n));
    modelX.col(4) = modelX.col(0)+modelX.col(1)+modelX.col(2)+modelX.col(5)+modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    modelX.col(7) = modelX.col(5)-modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    // Rcpp::Rcout << "Model defined  "<< modelX<<std::endl; 
    // Rcpp::Rcout << "Model defined  "<< resY<<std::endl; 
    // Here we have the model defined for iteration #s
    
    ///////////////////////////////////////////////////////////////////////////
    //Then we start the for loop to run over iterations
    vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
    // X.row(0)=1;
    // X.row(1)=1;
    // X.row(2)=1;
    // Rcpp::Rcout << "Creates initial vector  "<< X<<std::endl;
    // Rcpp::Rcout << "Check modes initial  "<< check_modes(X)<<std::endl;
    int curr_temp = 0; // We start in the very first temperature (which should be 1)
    // Rcpp::Rcout << "Defines initial temperature "<< curr_temp<<std::endl;
    vec logpsi(t, fill::zeros); // Initialize vector logpsi, as many entries as temperatures
    //Starts simulation
    for(int i=0;i<numiter;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 
      
      
      double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
      
      int total_neighbors = p+t; // total number of neighbors is p spacial + t temperature
      // Rcpp::Rcout << "p= "<< p <<std::endl;
      // Rcpp::Rcout << "t= "<< t <<std::endl;
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(sum(current_coord)==0){ // If the current state is all zeros
        logpi_current=logL_0(resY);
      }else{ // If the current state is not all zeroes
        logpi_current = logLikelihood(modelX,resY,current_coord);
      }
      ////////////      
      //Compute weight for spatial neighbors
      double temporal=0;
      vec newX;
      for(int j=0; j<p;j++){
        // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
        newX = X;
        newX.row(j) = 1-X.row(j);
        //Rcpp::Rcout << newX << std::endl;
        uvec coord = find(newX==1);
        //Rcpp::Rcout << coord << std::endl;
        if(sum(coord)==0){// In case the state visited is all 0s
          temporal=(logL_0(resY)-logpi_current)*temperature;
        }else{// For every other state that is not all 0s
          temporal=((logLikelihood(modelX,resY,coord)-logpi_current)*temperature);
        }
//Apply balancing function to spatial neighbors /////////////////////////////
        if(curr_temp<J){//For all other temps
          probs(j)=(temporal/2) - log(p_double);
        }else{//for temperature J
          if(temporal<0){probs(j) = temporal - log(p_double);}else{probs(j) =-log(p_double);}   
        }
        
      }
      ////////////      
      // Compute weight for temperature neighbors
      double temp_nei = 0;
      for(int j=p; j<(p+t);j++){
        // Rcpp::Rcout << "Checks temp neighbors  "<< j<<std::endl; 
        // Rcpp::Rcout << "Checks next neighbor  "<< probs(j+1)<<std::endl; 
        // Rcpp::Rcout << "index for temp  "<< j-p<<std::endl; 
        temp_nei = conv_to<double>::from(temp.row(j-p)); //Get the value of the temperature
        // Rcpp::Rcout << "temperature  "<< temp.row(j-p)<<std::endl;
        //Usually we deleted the entry for the current temperature
        //But it's easier just to assign a probability 0 to choose it
        if(temp_nei==temperature){probs(j)=0;}else{//
          // temporal=(logpi_current*(temp_nei-temperature)); //This was without the logpsi factors
          temporal=(logpi_current*(temp_nei-temperature) + conv_to<double>::from(logpsi.row(j-p))-conv_to<double>::from(logpsi.row(curr_temp)));
          //Apply balancing function to temperature neighbors /////////////////////////////           
          if(temporal<0){probs(j) = temporal - log(J);}else{probs(j)=-log(J);}
        }
        // Rcpp::Rcout << "Finish iteration of loop  "<< j<<std::endl; 
      }
      
      //Choose the next neighbor
      vec u = Rcpp::runif(p+t);
      vec probs_choose = log(-log(u)) - probs;
      
      //Find the index of the minimum element. source:https://gallery.rcpp.org/articles/vector-minimum/
      //This corresponds to choosing that neighbor
      int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
      // Rcpp::Rcout << neigh_pos << std::endl;
      
      if(neigh_pos<p){
        X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
      }else{
        curr_temp = neigh_pos-p; //Update the temperature
      }
      
      // Update logpsi vector
      logpsi = update_logpsi(logpsi,curr_temp,double(i),t-1);
      // Rcpp::Rcout << logpsi << std::endl;
      //update_logpsi(vec logpsi, int curr_temp, double iteration, int J)
      
      // Store the new values to report
      // modes_visited.row((s*numiter)+i)=check_modes(X)+1;
      modes_visited.submat((s*numiter)+i,0,(s*numiter)+i,0)=check_modes(X)+1;
      modes_visited.submat((s*numiter)+i,1,(s*numiter)+i,1)=curr_temp+1;
    }// End of for loop for iterations
    // Rcpp::Rcout << "Finish simulation "<< s << std::endl;
    // Rcpp::Rcout << "X= "<< X << std::endl;
  }//End of for loop for simulations
  return modes_visited;
}//End of function

// This function reports which modes were visited in which iteration and the temperature
// This is for model 3 which uses min only in temperature J and sq balancing function for every other
// [[Rcpp::export]]
vec Simulation_mod_IIT(int n, int p, int numsim, int numiter){
  // Attempt to create a single function that does everything
  
  //////////////////////////////////////////////////
  vec modes_visited(numiter * numsim);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  for(int s=0; s<numsim;s++){
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
    //// Creates model
    // Create the model matrix X and response Y
    mat modelX(n,p); // n rows and p columns (p corresponds to the dimension of the problem)
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
    
    //Create multi modality
    modelX.col(3) = modelX.col(1)-modelX.col(2)+as<vec>(Rcpp::rnorm(n));
    modelX.col(4) = modelX.col(0)+modelX.col(1)+modelX.col(2)+modelX.col(5)+modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    modelX.col(7) = modelX.col(5)-modelX.col(6)+as<vec>(Rcpp::rnorm(n));
    // Rcpp::Rcout << "Model defined  "<< modelX<<std::endl; 
    // Rcpp::Rcout << "Model defined  "<< resY<<std::endl; 
    // Here we have the model defined for iteration #s
    
    ///////////////////////////////////////////////////////////////////////////
    //Then we start the for loop to run over iterations
    vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
    //Starts simulation
    for(int i=0;i<numiter;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s << " Iteration: " << i << std::endl;}
      int total_neighbors = p; // total number of neighbors is p spacial + t temperature
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(sum(current_coord)==0){ // If the current state is all zeros
        logpi_current=logL_0(resY);
      }else{ // If the current state is not all zeroes
        logpi_current = logLikelihood(modelX,resY,current_coord);
      }
      ////////////      
      //Compute weight for spatial neighbors
      double temporal=0;
      vec newX;
      for(int j=0; j<p;j++){
        // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
        newX = X;
        newX.row(j) = 1-X.row(j);
        //Rcpp::Rcout << newX << std::endl;
        uvec coord = find(newX==1);
        //Rcpp::Rcout << coord << std::endl;
        if(sum(coord)==0){// In case the state visited is all 0s
          temporal=(logL_0(resY)-logpi_current);
        }else{// For every other state that is not all 0s
          temporal=((logLikelihood(modelX,resY,coord)-logpi_current));
        }
        //Apply balancing function to spatial neighbors /////////////////////////////
          probs(j)=exp(temporal/2)/p_double;
      }
      ////////////      
 
      //Choose the next neighbor
      vec u = Rcpp::runif(p);
      vec probs_choose = -log(u)/probs;
      
      //Find the index of the minimum element. source:https://gallery.rcpp.org/articles/vector-minimum/
      //This corresponds to choosing that neighbor
      int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
      // Rcpp::Rcout << neigh_pos << std::endl;
      
      X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor

      // Store the new values to report
      modes_visited.row((s*numiter)+i)=check_modes(X)+1;

    }// End of for loop for iterations
    // Rcpp::Rcout << "Finish simulation "<< s << std::endl;
    // Rcpp::Rcout << "X= "<< X << std::endl;
  }//End of for loop for simulations
  return modes_visited;
}//End of function