#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec fast_conv(const arma::vec& x, const arma::vec& y) {
  return arma::conv(x, y);
}

// [[Rcpp::export]]
Rcpp::IntegerVector detect_peaks(const arma::vec& signal, 
                       const double lower_bound, 
                       const int refractory) {
  Rcpp::IntegerVector out;
  int out_n = 0;
  
  for (int i = 1; i < signal.n_elem; i++) 
  {
    if (out_n > 0 && (i - out[out_n - 1]) < refractory) continue;
    
    if (signal[i] > lower_bound && 
        signal[i] > signal[i - 1] && 
        signal[i] > signal[i + 1])
    {
      out.push_back(i);
      out_n += 1;
    }
  }
  return out + 1;
}