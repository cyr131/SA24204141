#include <Rcpp.h>
using namespace Rcpp;

//' @title EX correlation matrix function
//' @description generate exchangable working correlation matrix
//' @param t Specify matrix dimensions
//' @param alpha Specify the correlation parameter
//' @return EX correlation matrix with t_dimensions and alpha_parameter
//' @examples
//' \dontrun{
//' t <- 5
//' alpha <- 0.5
//' exchCpp(t,alpha)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix exchCpp(int t, double alpha) {
  NumericMatrix exch(t, t);

  // 初始化对角线元素为1
  for (int i = 0; i < t; ++i) {
    exch(i, i) = 1.0;
  }

  // 设置非对角线元素为alpha
  for (int i = 0; i < t; ++i) {
    for (int j = 0; j < t; ++j) {
      if (i != j) {
        exch(i, j) = alpha;
      }
    }
  }

  return exch;
}


//' @title AR1 correlation matrix function
//' @description generate autoregression working correlation matrix
//' @param t Specify matrix dimensions
//' @param alpha Specify the correlation parameter
//' @return AR1 correlation matrix with t_dimensions and alpha_parameter
//' @examples
//' \dontrun{
//' t <- 5
//' alpha <- 0.5
//' ar1Cpp(t,alpha)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix ar1Cpp(int t, double alpha) {
   NumericMatrix ar1(t, t);
   for(int i = 0; i < t; ++i) {
     for(int j = i; j < t; ++j) {
       if(i == j) {
         ar1(i, j) = 1;
       } else {
         ar1(i, j) = pow(alpha, j - i);
         ar1(j, i) = ar1(i, j);
       }
     }
   }
   return ar1;
 }

//' @title ST correlation matrix function
//' @description generate stationary working correlation matrix
//' @param t Specify matrix dimensions
//' @param a Specify non-diagonal element
//' @return ST correlation matrix with t_dimensions and non-diagonal element a
//' @examples
//' \dontrun{
//' t <- 6
//' a <- c(0.5,0.4,0.3,-0.2,-0.1)
//' statCpp(t,a)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix statCpp(int t, NumericVector a) {
  int m = a.size();
  NumericMatrix stat(t, t);

  // 设置对角线元素为1
  for(int i = 0; i < t; ++i) {
    stat(i, i) = 1;
  }

  // 设置非对角线元素
  for(int i = 0; i < t; ++i) {
    for(int j = 0; j < t; ++j) {
      if(abs(i - j) > 0 && abs(i - j) <= m) {
        stat(i, j) = a[abs(i - j) - 1];
      }
    }
  }

  return stat;
}

