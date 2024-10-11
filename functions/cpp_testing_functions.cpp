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
mat test_mat(mat M, vec X, int p){
  for(int i=0;i<p;i++){
    X.row(i)=1-X.row(i);
    M.row(i)=X.t();
  }
  return M;
}
// [[Rcpp::export]]
vec test_vec(int t){
  vec vector(t,fill::zeros);
  
  for(int i=0;i<t;i++){
    vector.row(i)=i+3;
  }
  return vector;
}

// [[Rcpp::export]]
uvec check_ones(vec X){
  uvec pos=find(X==1);
  return pos;
}


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



// Function to read a matrix from a file providing the file's name
// [[Rcpp::export]]
mat read_file_cpp(String filename) {
  //std::stringstream filename
  // Rcpp::Rcout << "Entra funcion" <<std::endl;
  mat modelX;
  
  bool status = modelX.load(filename, arma::csv_ascii);
  // Rcpp::Rcout << "status= "<< status <<std::endl;
  // Rcpp::Rcout << "X= "<< modelX <<std::endl;
  if (!status) {
    String s("Error loading ");
    s+=filename;
    Rcpp::stop(s);
  }
  return modelX;
}

// [[Rcpp::export]]
vec read_Y_cpp(String filename) {
  //std::stringstream filename
  // Rcpp::Rcout << "Entra funcion" <<std::endl;
  vec resY;
  
  bool status = resY.load(filename, arma::csv_ascii);
  // Rcpp::Rcout << "status= "<< status <<std::endl;
  // Rcpp::Rcout << "X= "<< modelX <<std::endl;
  if (!status) {
    String s("Error loading ");
    s+=filename;
    Rcpp::stop(s);
  }
  return resY;
}

// [[Rcpp::export]]
vec read_Y_i(int i) {
  //std::stringstream filename
  // Rcpp::Rcout << "Entra funcion" <<std::endl;
  vec resY;
  std::string filename = "models/resY" + std::to_string(i) + ".csv";
  bool status = resY.load(filename, arma::csv_ascii);
  // Rcpp::Rcout << "status= "<< status <<std::endl;
  // Rcpp::Rcout << "X= "<< modelX <<std::endl;
  if (!status) {
    Rcpp::stop("Error loading file: " + filename);
  }
  return resY;
}

// [[Rcpp::export]]
double read_file_inside_function(int n){
  
  vec X={1,1,0,1};
  uvec coord = find(X==1);
  Rcpp::Rcout << coord <<std::endl;
  mat modelX;
  vec resY;
  double ll;
  bool status;
    // String filename_model("models/modelX");
  // String filename_res("models/resY");
  // String temp("");
  // String contador;
  // std::string filename_model="C:/Users/ralex/Documents/src/comp_report/codes/modelX";
  // std::string filename_res="C:/Users/ralex/Documents/src/comp_report/codes/resY";
  // std::string temp;
  for(int i=0;i<n;i++){
    Rcpp::Rcout << "For loop iter= "<< i <<std::endl;
    // contador= std::to_string(i+1) + ".csv";
    // temp= filename_model += contador;
    // Rcpp::Rcout << temp.get_cstring() <<std::endl;
    status = modelX.load("models/modelX" + std::to_string(i+1) + ".csv", arma::csv_ascii);
     Rcpp::Rcout << status <<std::endl;
    // Rcpp::Rcout << modelX <<std::endl;
    // temp= filename_res += contador;
    // Rcpp::Rcout << temp <<std::endl;
    status = resY.load("models/resY" + std::to_string(i+1) + ".csv", arma::csv_ascii);
    Rcpp::Rcout << status <<std::endl;
    ll=logLikelihood(modelX,resY,coord);
    Rcpp::Rcout << ll <<std::endl;
  }
  return ll;
}




//////////////////////////////////
//Some useful links

// How to use functions inside other functions
//If they're in the same file, they can interact
//https://stackoverflow.com/questions/45296381/can-we-use-rcpp-with-multiple-c-functions/45296773#45296773


//Translating code from R to armadillo



// COnvert vectors from Rcpp to armadillo
// https://statisticsglobe.com/convert-vector-from-rcpparmadillo-to-rcpp-r

// To subset from a vector, 
// If i'm using uvec I need to use v.elem()
//If im using int I need to use v.row()

