#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;


// Function to orthogonalize vector 'x' with respect to matrix 'A'
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd orthogonalization (Eigen::VectorXd x, const Eigen::MatrixXd& A) {
  int ncol_A = A.cols();

  for (int i = 0; i < ncol_A; ++i) {
    double proj_coeff = x.dot(A.col(i)) / A.col(i).squaredNorm();
    x -= proj_coeff * A.col(i);
  }

  return x;
}


// Function to thick restart: Truncate and then expand to Hm
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List truncate_and_expand_Cpp (const Eigen::MatrixXd& A, Eigen::MatrixXd V, Eigen::VectorXd u,
                         Eigen::VectorXd b1, Eigen::MatrixXd Tm, Eigen::MatrixXd Sm,
                         int n, int m, int k, double tol0 = 1e-8) {
  // eigenvalues are already reordered by SchurDecom()
  // Truncate to a Krylov-Schur decomposition of order k
  Eigen::MatrixXd V11 = (V * Sm).leftCols(k);
  Eigen::MatrixXd T11 = Tm.topLeftCorner(k, k);
  Eigen::VectorXd u11 = u / u.norm();
  Eigen::VectorXd b11 = (Sm.transpose() * b1).topRows(k) * u.norm();
  b11 = (b11.array().abs() >= tol0 * Eigen::ArrayXd::Ones(b11.size())).select(b11, 0.0);

  // Extend to a Krylov decomposition of order m
  Eigen::MatrixXd Unew(n, m);
  Unew.block(0, 0, n, k + 1) << V11, u11;

  Eigen::MatrixXd Snew = Eigen::MatrixXd::Zero(m, m);
  Snew.block(0, 0, k, k) << T11;
  Snew.block(k, 0, 1, k) << b11.transpose();
  Snew.block(0, k, k, 1) << b11;

  Eigen::VectorXd z = A * u11;
  double alpha = u11.dot(z);
  Eigen::VectorXd rr = z - alpha * u11 - V11 * b11;
  if (rr.norm() < std::sqrt(alpha * alpha + b11.squaredNorm())) {
    rr = orthogonalization(rr, Unew.block(0, 0, n, k + 1));
  }

  double beta = rr.norm();
  Snew(k, k) = alpha;

  for (int j = k + 1; j < m; j++) {
    Snew(j - 1, j) = Snew(j, j - 1) = beta;
    Unew.col(j) = rr / beta;
    z = A * Unew.col(j);
    alpha = Unew.col(j).dot(z);
    rr = z - alpha * Unew.col(j) - beta * Unew.col(j - 1);
    if (rr.norm() < std::sqrt(alpha * alpha + beta * beta)) {
      rr = orthogonalization(rr, Unew.block(0, 0, n, j + 1));
    }
    beta = rr.norm();
    Snew(j, j) = alpha;
  }

  Eigen::MatrixXd u_new = rr / beta;
  Eigen::MatrixXd b_new = (Eigen::VectorXd(m) << Eigen::VectorXd::Zero(m - 1), beta).finished();
  return List::create(Named("U") = Unew, Named("S") = Snew,
                      Named("u") = u_new, Named("b") = b_new.transpose());
}


// Function to perform Householder transformation on a matrix
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List HouseholderTransform_Cpp (int start, Eigen::MatrixXd& S1, Eigen::MatrixXd& V) {
  int m = S1.rows();

  for (int j = start; j < m-2; j++) {
    Eigen::VectorXd v = S1.block(j + 1, j, m - j, 1);

    if ((v.array() == 0).all()) {
      continue;
    }

    double sigma = S1.coeff(j + 1, j) > 0 ? v.norm() : -v.norm();
    Eigen::VectorXd unitVector = Eigen::VectorXd::Zero(m - j);
    unitVector[0] = 1.0;
    v += sigma * unitVector;
    Eigen::MatrixXd V1 = Eigen::MatrixXd::Identity(m, m);
    V1.block(j + 1, j + 1, m - j - 1, m - j - 1) -= 2 * v * v.transpose() / v.squaredNorm();
    S1 = V1.transpose() * S1 * V1;
    S1.block(j, std::min(j + 2, m), 1, m - std::min(j + 2, m)).setZero();
    S1.block(std::min(j + 2, m), j, m - std::min(j + 2, m), 1).setZero();
    V *= V1;
  }

  return List::create(Named("S1") = S1, Named("V") = V);
}
