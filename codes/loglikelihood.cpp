#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
vec logLikelihood_num(vec res) {
  uword n = res.n_elem;
  
  vec rss=res.t()*res;
  vec logl = -0.5*(n*log(2*M_PI) + n*log(rss) -n*log(n) + n);
  return logl;
}



/*** R
# logLikelihood(c(1,2,3))
# logLikelihood(c(1,2,3))
*/
