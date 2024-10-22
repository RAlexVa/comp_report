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
//Adding a penalization for the number of terms considering a normal prior for the betas
//add a penalization for small coefficients and number of parameters
logl=logl - 0.5 *arma::dot(coef,coef) - 3*log(n)*(subX.n_cols+1) + sum(log(abs(coef)));

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
mat Simulation_mod1(int n, int p, int startsim,int endsim, int numiter, vec temp,int t){
  // This function takes the # of the starting and ending simulation so we can run a specific number of simulations

  //////////////////////////////////////////////////
  int total_sim = (endsim-startsim+1);
  mat modes_visited(numiter *total_sim ,2);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  double t_double = double(t);//
  // Rcpp::Rcout << "t= "<< t <<std::endl;
  // Rcpp::Rcout << "t_double= "<< t_double <<std::endl;
  double J=t_double-1;//Number of temperatures minus 1
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  // s always has to start from 0 since at the end of the loop we use the varaible s to store values in a matrix
  for(int s=0; s<total_sim;s++){ //The loop considers the start and ending simulation #s
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
//// Creates model
    // Read the model matrix X and response Y according to the # of simulation
    mat modelX;
    vec resY;
    bool status; //To check if there are any issues with the reading
    status = modelX.load("models/modelX" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    // Rcpp::Rcout << "Reads file  "<< s+startsim<<std::endl;
    if (!status) {Rcpp::stop("Error loading file: ModelX for simulation" + std::to_string(s+startsim));}
    status = resY.load("models/resY" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    if (!status) {Rcpp::stop("Error loading file: ResY for simulation" + std::to_string(s+startsim));}
    // Here we have the model defined for iteration #s
    // vec Xtest={1,1,0,1};
    // uvec coord = find(Xtest==1);
    // Rcpp::Rcout << "Likelihood model  "<< s+startsim<< " is "<< logLikelihood(modelX,resY,coord)<<std::endl;
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
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 


      double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
      
      int total_neighbors = p+t; // total number of neighbors is p spacial + t temperature
      // Rcpp::Rcout << "p= "<< p <<std::endl;
      // Rcpp::Rcout << "t= "<< t <<std::endl;
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(current_coord.empty()){ // If the current state is all zeros
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
        if(coord.empty()){// In case the state visited is all 0s
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
mat Simulation_mod2(int n, int p, int startsim,int endsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything
  
  //////////////////////////////////////////////////
  int total_sim = (endsim-startsim+1);
  mat modes_visited(numiter * total_sim,2);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  double t_double = double(t);//
  // Rcpp::Rcout << "t= "<< t <<std::endl;
  // Rcpp::Rcout << "t_double= "<< t_double <<std::endl;
  double J=t_double-1;//Number of temperatures minus 1
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  for(int s=0; s<total_sim;s++){
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
//// Creates model
    // Read the model matrix X and response Y according to the # of simulation
    mat modelX;
    vec resY;
    bool status; //To check if there are any issues with the reading
    status = modelX.load("models/modelX" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    // Rcpp::Rcout << "Reads file  "<< s+startsim<<std::endl;
    if (!status) {Rcpp::stop("Error loading file: ModelX for simulation" + std::to_string(s+startsim));}
    status = resY.load("models/resY" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    if (!status) {Rcpp::stop("Error loading file: ResY for simulation" + std::to_string(s+startsim));}
    // Here we have the model defined for iteration #s   
    // vec Xtest={1,1,0,1};
    // uvec coord = find(Xtest==1);
    // Rcpp::Rcout << "Likelihood model  "<< s+startsim<< " is "<< logLikelihood(modelX,resY,coord)<<std::endl;
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
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 
      
      
      double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
      
      int total_neighbors = p+t; // total number of neighbors is p spacial + t temperature
      // Rcpp::Rcout << "p= "<< p <<std::endl;
      // Rcpp::Rcout << "t= "<< t <<std::endl;
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(current_coord.empty()){ // If the current state is all zeros
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
        if(coord.empty()){// In case the state visited is all 0s
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
mat Simulation_mod3(int n, int p, int startsim, int endsim, int numiter, vec temp,int t){
  // Attempt to create a single function that does everything
  
  //////////////////////////////////////////////////
  int total_sim=(endsim-startsim+1);
  mat modes_visited(numiter * total_sim,2);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  double t_double = double(t);//
  // Rcpp::Rcout << "t= "<< t <<std::endl;
  // Rcpp::Rcout << "t_double= "<< t_double <<std::endl;
  double J=t_double-1;//Number of temperatures minus 1
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  for(int s=0; s<total_sim;s++){
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
  //// Creates model
    // Read the model matrix X and response Y according to the # of simulation
    mat modelX;
    vec resY;
    bool status; //To check if there are any issues with the reading
    status = modelX.load("models/modelX" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    // Rcpp::Rcout << "Reads file  "<< s+startsim<<std::endl;
    if (!status) {Rcpp::stop("Error loading file: ModelX for simulation" + std::to_string(s+startsim));}
    status = resY.load("models/resY" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    if (!status) {Rcpp::stop("Error loading file: ResY for simulation" + std::to_string(s+startsim));}
    // Here we have the model defined for iteration #s   
    // vec Xtest={1,1,0,1};
    // uvec coord = find(Xtest==1);
    // Rcpp::Rcout << "Likelihood model  "<< s+startsim<< " is "<< logLikelihood(modelX,resY,coord)<<std::endl;
    
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
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
      // Rcpp::Rcout << "Starts iteration  "<< i<<std::endl; 
      
      
      double temperature = conv_to<double>::from(temp.row(curr_temp)); // Current temperature
      
      int total_neighbors = p+t; // total number of neighbors is p spacial + t temperature
      // Rcpp::Rcout << "p= "<< p <<std::endl;
      // Rcpp::Rcout << "t= "<< t <<std::endl;
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(current_coord.empty()){ // If the current state is all zeros
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
        if(coord.empty()){// In case the state visited is all 0s
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
vec Simulation_mod_IIT(int n, int p, int startsim, int endsim, int numiter){
  // Attempt to create a single function that does everything
  
  //////////////////////////////////////////////////
  int total_sim = (endsim-startsim+1);
  vec modes_visited(numiter * total_sim);//Matrix to store the modes visited and temperature
  double p_double = double(p);// 
  // Repeats according to the number of simulations
  /////////////////////////////////////////////////////
  for(int s=0; s<total_sim;s++){
    //Rcpp::Rcout << "Starts simulation number  "<< s<<std::endl;
  //// Creates model
    // Read the model matrix X and response Y according to the # of simulation
    mat modelX;
    vec resY;
    bool status; //To check if there are any issues with the reading
    status = modelX.load("models/modelX" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    // Rcpp::Rcout << "Reads file  "<< s+startsim<<std::endl;
    if (!status) {Rcpp::stop("Error loading file: ModelX for simulation" + std::to_string(s+startsim));}
    status = resY.load("models/resY" + std::to_string(s+startsim) + ".csv", arma::csv_ascii);
    if (!status) {Rcpp::stop("Error loading file: ResY for simulation" + std::to_string(s+startsim));}
    // Here we have the model defined for iteration #s   
    // vec Xtest={1,1,0,1};
    // uvec coord = find(Xtest==1);
    // Rcpp::Rcout << "Likelihood model  "<< s+startsim<< " is "<< logLikelihood(modelX,resY,coord)<<std::endl;
    
    ///////////////////////////////////////////////////////////////////////////
    //Then we start the for loop to run over iterations
    vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
    //Starts simulation
    for(int i=0;i<numiter;i++){
      if (i % 1000 == 1) {Rcpp::Rcout << "Simulation: " << s+startsim << " Iteration: " << i << std::endl;}
      int total_neighbors = p; // total number of neighbors is p spacial + t temperature
      vec probs(total_neighbors, fill::zeros); //probabilities
      // Compute likelihood of the current state
      double logpi_current=0;
      uvec current_coord = find(X==1);
      if(current_coord.empty()){ // If the current state is all zeros
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
        if(coord.empty()){// In case the state visited is all 0s
          temporal=(logL_0(resY)-logpi_current);
        }else{// For every other state that is not all 0s
          temporal=((logLikelihood(modelX,resY,coord)-logpi_current));
        }
        //Apply balancing function to spatial neighbors /////////////////////////////
          probs(j)=(temporal/2) - log(p_double);
      }
      ////////////      
 
      //Choose the next neighbor
      vec u = Rcpp::runif(p);
      vec probs_choose = log(-log(u)) - probs;
      
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

// [[Rcpp::export]]
mat readmodelX(int sim){
  mat modelX;
  bool status; //To check if there are any issues with the reading
  status = modelX.load("models/modelX" + std::to_string(sim) + ".csv", arma::csv_ascii);
  // Rcpp::Rcout << "Reads file  "<< s+startsim<<std::endl;
  if (!status) {Rcpp::stop("Error loading file: ModelX for simulation" + std::to_string(sim));}
  return modelX;
}

// [[Rcpp::export]]
vec readY(int sim){
  vec resY;
  bool status; //To check if there are any issues with the reading
  status = resY.load("models/resY" + std::to_string(sim) + ".csv", arma::csv_ascii);
  if (!status) {Rcpp::stop("Error loading file: ResY for simulation" + std::to_string(sim));}
  return resY;
}

// [[Rcpp::export]]
double bf_sq(double x){
  return x/2;
}

// [[Rcpp::export]]
double bf_min(double x){
  double result;
  if(x<0){
    result = x;
  }else{
    result = 0;
  }
  return result;
}

//This function we don't export to R because it generates errors
double invoke(double x, double (*func)(double)) {
  return func(x);
}

// [[Rcpp::export]]
double bal_func(double x,String chosen){
  if (chosen == "sq") {
    return invoke(x, &bf_sq);
  } else if (chosen == "min") {
    return invoke(x, &bf_min);
  } else {
    cout << "Unknown operation!" << endl;
    Rcpp::Rcout <<"Unknown operation!" << std::endl;
    return 0; // Default return for unknown operation
  }
}

////////////////////////////////////////////////////

// [[Rcpp::export]]
vec RF_update(vec X, String chosen_bf, mat modelX, vec resY){
  int total_neighbors = X.n_rows; // total number of neighbors is p spacial
  vec probs(total_neighbors, fill::zeros); //probabilities
  // Compute likelihood of the current state
  double logpi_current=0;
  uvec current_coord = find(X==1);
  if(current_coord.empty()){ // If the current state is all zeros
    logpi_current=logL_0(resY);
  }else{ // If the current state is not all zeroes
    logpi_current = logLikelihood(modelX,resY,current_coord);
  }
  Rcpp::Rcout << "current loglik: "<< logpi_current << std::endl;
  ////////////      
  //Compute weight for all neighbors
  double temporal=0;
  vec newX;
  for(int j=0; j<total_neighbors;j++){
    // Rcpp::Rcout << "Starts checking neighbors  "<< j<<std::endl; 
    newX = X;
    newX.row(j) = 1-X.row(j);
    //Rcpp::Rcout << newX << std::endl;
    uvec coord = find(newX==1);
    //Rcpp::Rcout << coord << std::endl;
    if(coord.empty()){// In case the state visited is all 0s
      temporal=logL_0(resY)-logpi_current;
    }else{// For every other state that is not all 0s
      temporal=logLikelihood(modelX,resY,coord)-logpi_current;
    }
    // Rcpp::Rcout << "loglik antes de bf: "<< temporal <<" en neighbor "<<j+1<< std::endl;
    //Apply balancing function to log probabilities /////////////////////////////
    probs(j)=bal_func(temporal, chosen_bf);
  }
  ////////////      
  
  //Choose the next neighbor
  vec u = Rcpp::runif(total_neighbors);
  vec probs_choose = log(-log(u)) - probs;
  
  //Find the index of the minimum element. 
  //This corresponds to choosing that neighbor
  int neigh_pos = (std::min_element(probs_choose.begin(), probs_choose.end()))-probs_choose.begin();
  // Rcpp::Rcout <<"probs vector: "<< probs << std::endl;
  // Rcpp::Rcout <<"chosen neighbor: "<< neigh_pos << std::endl;
  
  X.row(neigh_pos) = 1-X.row(neigh_pos); //modify the coordinate of the chosen neighbor
  return X;
}


// [[Rcpp::export]]
vec RF_IIT_sim(int numsim, int numiter,int p){
  vec X(p,fill::zeros); // The starting state of all simulations is a vector full of zeroes
  
  for(int s=0;s<numsim;s++){
    mat modelX=readmodelX(s+1);
    vec resY=readY(s+1);
    for(int i=0;i<numiter;i++){
      X = RF_update(X,"sq",modelX,resY);
    }
    Rcpp::Rcout <<"Final state "<< X << std::endl;
  }
  return X;
}


