#' A implementation of Arnoldi Decomposition
#' @param A input matrix, to decompose it into AV = VH + r*t(e)
#' @param b vector for Krylov-subspace
#' @param V0 start V
#' @param H0 start H
#' @param r0 start r
#' @param n dimension of A
#' @param m dimension of H
#' @param k need k largest eigenvalues
#' @param First whether it is the first Arnoldi, before any restart
#' @return A list containing H, V and r
#' @import Rcpp
#' @export
#'
Arnoldi <- function(A, b, V0, H0, r0, n, m, k, First=FALSE) {
  # Hessenberg Decomposition AV = VH, where VV' = I, H is upper Hessenberg
  V <- matrix(0, nrow = n, ncol = m + 1)
  H <- matrix(0, nrow = m + 1, ncol = m)
  if (First) {
    V[, 1] <- b / norm2(b)
    start <- 2
  } else {
    H[1:k, 1:k] <- H0[1:k, 1:k]
    V[, 1:k] <- V0[, 1:k]
    r0_orth <- r0
    for (l in 1:k) {
      proj <- r0 %*% V[, l]
      r0_orth <- r0_orth - proj[1,1] * as.vector(V[, l])
    }
    H[k+1, k] <- norm2(r0_orth)
    V[, k+1] <- r0_orth / H[k+1, k]
    start <- k + 2
  }

  # Arnoldi Algorithm
  ArnoldiOut <- Arnoldi_Cpp(A, V, H, start, n, m)
  # A %*% ArnoldiOut$V - ArnoldiOut$V %*% ArnoldiOut$H - ArnoldiOut$r %*% t(c(rep(0,m-1),1))
  return(list(H=ArnoldiOut$H, V=ArnoldiOut$V, r=ArnoldiOut$r))
}


#' A implementation of checking convergence of Arnoldi Decomposition
#' @param A input matrix, to decompose it into AV = VH + r*t(e)
#' @param b vector for Krylov-subspace
#' @param V0 start V
#' @param H0 start H
#' @param r0 start r
#' @param n dimension of A
#' @param m dimension of H
#' @param k need k largest eigenvalues
#' @param First whether it is the first Arnoldi, before any restart
#' @param tol tolerance of convergence
#' @return A list containing convergence, H, V, r or eigenpairs
#'
Arnoldi_converge <- function(A, b, V0, H0, r0, n, m, k, First=FALSE, tol=1e-6) {
  ArnoldiVH <- Arnoldi(A, b, V0, H0, r0, n, m, k, First)
  H <- ArnoldiVH$H
  V <- ArnoldiVH$V
  r <- ArnoldiVH$r

  # Solve the eigenpairs of Hessenberg Matrix H
  SchurH <- SchurDecom(H)
  # Ritz Vector conditions for termination
  Ritz <- A %*% V %*% SchurH$Q - t(diag(SchurH$T) * t(V %*% SchurH$Q))
  if (max(apply(Ritz[, 1:k], 2, norm2)) < tol*n) {
    return(list(convergence=TRUE, eigenVal=diag(SchurH$T), eigenVec=V %*% SchurH$Q))
  } else {
    return(list(convergence=FALSE, eigenVal=diag(SchurH$T), H=H, V=V, r=r, r=r))
  }
}


#' A implementation of Implicitly Restart Arnoldi
#' @param H Hessenberg matrix, H of AV = VH + r*t(e)
#' @param V orthogonal matrix, V of AV = VH + r*t(e)
#' @param r vector, r of AV = VH + r*t(e)
#' @param eigenVal eigenvalues of H
#' @param k need k largest eigenvalues
#' @param m dimension of H
#' @return A list containing H, V and r after restarting
#' @import Rcpp
#'
ImplicitlyRestart <- function(H, V, r, eigenVal, k, m) {
  H1 <- H
  V1 <- V
  G <- diag(m)

  # restart by removing unwanted eigenvalues
  for (i in (k+1):m) {
    sigma <- eigenVal[i]
    qrH <- GivensRotation_Cpp(H1, sigma, m)
    H1 <- qrH$H_new
    for (i in 3:m) { H1[i, 1:(i-2)] <- 0 }
    for (i in 1:(m-2)) { H1[i, (i+2):m] <- 0 }
    V1 <- V1 %*% t(qrH$G)
    G <- G %*% t(qrH$G)
  }
  r1 <- V1[, k+1] * H1[k+1, k] + r * G[m, k]
  # A %*% V1[, 1:k] - V1[,1:k] %*% H1[1:k, 1:k] - r1 %*% t(c(rep(0, k-1),1))

  return(list(H=H1[1:k, 1:k], V=V1[, 1:k], r=r1))
}


#' A implementation IRAM for EVD
#' @param A input matrix
#' @param k need k largest eigenvalues
#' @param max_iter maximum of iteration number
#' @return A list containing eigenvalues and eigenvectors
#' @export
#'
eigen_IRAM <- function(A, k, max_iter=10000) {
  # check whether A is a square matrix
  if (nrow(A) != ncol(A)) { stop('Input must be a matrix!') }
  n <- nrow(A)
  m <- min(2 * k, n)
  # randomly initialize b
  b <- runif(n, -1, 1)

  eigenH <- Arnoldi_converge(A, b, V0=NULL, H0=NULL, r0=NULL, n, m, k, First=TRUE)

  iter <- 0
  while (iter < max_iter) {
    if (eigenH$convergence == TRUE) {
      return(list(eigenVal=eigenH$eigenVal[1:k], eigenVec=eigenH$eigenVec[, 1:k], iter=iter))
    } else {
      iter <- iter + 1
      newVHr <- ImplicitlyRestart(eigenH$H, eigenH$V, eigenH$r, eigenH$eigenVal, k, m)
      eigenH <- Arnoldi_converge(A, b, newVHr$V, newVHr$H, newVHr$r, n, m, k, First=FALSE)
    }
  }
  stop('IRAM not convergence!')
}
