#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;


// Function to complete Arnoldi Rotation in Cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List Arnoldi_Cpp (const Eigen::MatrixXd& A, Eigen::MatrixXd V, Eigen::MatrixXd H, int start, int n, int m) {
  Eigen::VectorXd r(n);

  for (int j = start-1; j < m+1; j++) {
    Eigen::VectorXd Vj = A * V.col(j - 1);
    for (int i = 0; i <= j-1; i++) {
      H(i, j - 1) = V.col(i).dot(Vj);
      Vj -= H(i, j - 1) * V.col(i);
    }
    H(j, j - 1) = Vj.norm();
    V.col(j) = Vj / H(j, j - 1);
  }

  r = V.col(m) * H(m, m-1);
  V = V.leftCols(m).eval();
  H = H.topRows(m).eval();

  // Use mask matrix to set elements outside the tridiagonal to zero
  Eigen::MatrixXd mask = Eigen::MatrixXd::Ones(n, m);
  mask.diagonal(0).setOnes();
  mask.diagonal(1).setOnes();
  mask.diagonal(-1).setOnes();

  H.array() *= mask.array();

  return List::create(Named("V") = V, Named("H") = H, Named("r") = r);
}
