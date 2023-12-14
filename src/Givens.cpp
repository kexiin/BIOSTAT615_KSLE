#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;


// Function to perform Givens rotation using RcppEigen
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List GivensRotation_Cpp(Eigen::MatrixXd H1, double sigma, int m) {
  Eigen::MatrixXd H2 = H1 - sigma * Eigen::MatrixXd::Identity(m, m);
  Eigen::MatrixXd G = Eigen::MatrixXd::Identity(m, m);

  for (int i = 0; i < m - 1; i++) {
    double r = H2.block(i, i, 2, 1).norm();
    double c = H2(i, i) / r;
    double s = H2(i + 1, i) / r;

    Eigen::MatrixXd Gi = Eigen::MatrixXd::Identity(m, m);
    Gi.block(i, i, 2, 2) << c, s, -s, c;

    H2 = Gi * H2;
    G = Gi * G;
  }

  Eigen::MatrixXd H_new = H2 * G.transpose() + sigma * Eigen::MatrixXd::Identity(m, m);

  return List::create(Named("H_new") = H_new, Named("G") = G);
}
