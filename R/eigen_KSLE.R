#' A implementation of reducing the the Krylov form S to Hessenberg form by Householder matrix
#' @param U n*m orthogonal matrix, U of AU = US + ub
#' @param S m*m matrix, S of AU = US + ub
#' @param u n*1 matrix, u of AU = US + ub
#' @param b 1*m matrix, b of AU = US + ub
#' @param m dimension of S
#' @param k need k largest eigenvalues
#' @param tol0 tolerance for setting small numbers to 0
#' @return A list containing U, S, u and b, where S is tridiagonal matrix
#' @import Rcpp
#'
Arnoldi_Krylov <- function(U, S, u, b, m, k, tol0=1e-8) {
  start <- 1
  while ((length(which(S[(start+1):m, start] != 0)) == 0) & (start < k)) { start <- start + 1 }
  if (start == k) {
    return(list(U=U, S=S, u=u, b=b))
  }

  # the start column
  if (length(which(S[(start+1):m, start] != 0)) == 1) {
    V <- diag(m)
    V <- rbind(V[1:start, ], V[(start+2):(k+1), ], V[(start+1), ], V[(k+2):m,])
    S1 <- t(V) %*% S %*% V
  } else {
    non_zero_indice <- (which(S[(start+1):m, start] != 0))[1]
    em <- rep(0, m-start)
    em[non_zero_indice] <- 1
    v <- S[(start+1):m, start] - sign(S[non_zero_indice + 1, start]) * norm2(S[(start+1):m, start]) * em
    V <- diag(m)
    V[(start+1):m, (start+1):m] <- diag(m-start) - 2 * v %*% t(v) / norm2(v)^2
    S1 <- t(V) %*% S %*% V
    S1[start, min(start+2, m):m] <- 0
    S1[min(start+2, m):m, start] <- 0
  }

  # other columns
  Householder_to_Hess <- HouseholderTransform_Cpp(start, S1, V)

  U1 <- U %*% Householder_to_Hess$V
  b1 <- b %*% Householder_to_Hess$V
  # A %*% U1 - U1 %*% S1 - u %*% b1
  return(list(U=U1, S=Householder_to_Hess$S1, u=u, b=b1))
}


#' A implementation of checking convergence of Krylov-Schur Decomposition
#' @param A input matrix
#' @param V the orthogonal matrix of AV = VS + ub
#' @param H the diagonal matrix of Schur decomposition of S
#' @param G the orthogonal matrix of Schur decomposition of S
#' @param n dimension of A
#' @param k need k largest eigenvalues
#' @param tol tolerance of convergence
#' @return A list containing convergence and eigenpairs
#'
check_converge <- function(A, V, H, G, n, k, tol=1e-6) {
  od <- order(-diag(H))
  eigenVal <- diag(H)[od]
  eigenVec <- G[, od]

  Ritz <- A %*% V %*% eigenVec - t(eigenVal * t(V %*% eigenVec))
  # print(max(apply(Ritz[, 1:k], 2, norm2)))
  if (max(apply(Ritz[, 1:k], 2, norm2)) < tol*n) {
    return(list(convergence=TRUE, eigenVal=eigenVal, eigenVec=V %*% eigenVec))
  } else {
    return(list(convergence=FALSE, eigenVal=eigenVal, eigenVec=V %*% eigenVec))
  }
}


#' A implementation Krylov-Schur Algorithm for EVD
#' @param A input matrix
#' @param k need k largest eigenvalues
#' @param max_iter maximum of iteration number
#' @return A list containing eigenvalues and eigenvectors
#' @import Rcpp
#' @export
#'
eigen_KSLE <- function(A, k, max_iter=10000) {
  # check whether A is a square matrix
  if (nrow(A) != ncol(A)) { stop('Input must be a matrix!') }
  n <- nrow(A)
  m <- min(2 * k, n)
  # randomly initialize b
  bvec <- runif(n, -1, 1)

  # start with Arnoldi decomposition to get initial Krylov decomposition
  ArnoldiHV <- Arnoldi(A, bvec, V0=NULL, H0=NULL, r0=NULL, n, m, k, First=TRUE)
  SchurD <- SchurDecom(ArnoldiHV$H)
  em <- c(rep(0, m-1), 1)
  expandedUS <- truncate_and_expand_Cpp(A, ArnoldiHV$V, ArnoldiHV$r, em, SchurD$T, SchurD$Q, n, m, k)
  HessUS <- Arnoldi_Krylov(expandedUS$U, expandedUS$S, expandedUS$u, expandedUS$b, m, k)
  SchurD <- SchurDecom(HessUS$S)
  converge <- check_converge(A, HessUS$U, SchurD$T, SchurD$Q, n, k)

  iter <- 1
  while (iter < max_iter) {
    if (converge$convergence == TRUE) {
      return(list(eigenVal=converge$eigenVal[1:k], eigenVec=converge$eigenVec[, 1:k], iter=iter))
    } else {
      iter <- iter + 1
      expandedUS <- truncate_and_expand_Cpp(A, HessUS$U, HessUS$u, HessUS$b, SchurD$T, SchurD$Q, n, m, k)
      HessUS <- Arnoldi_Krylov(expandedUS$U, expandedUS$S, expandedUS$u, expandedUS$b, m, k)
      SchurD <- SchurDecom(HessUS$S)
      converge <- check_converge(A, HessUS$U, SchurD$T, SchurD$Q, n, k)
    }
  }
  stop('KSLE not convergence!')
}

