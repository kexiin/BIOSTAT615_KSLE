#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;


// Function to calculate norm2 for vectors
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double norm2 (const Eigen::VectorXd& x) {
  return x.norm();
}
