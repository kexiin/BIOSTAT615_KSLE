# library(control)

SchurDecom <- function(H, m) {
  qrH <- shiftedQR(H, m, check='strong')
  H1 <- qrH$H
  G1 <- qrH$G

  # H = G1 %*% H1 %*% t(G1)
  return(list(H=H1, G=G1))
}

truncate_and_expand <- function(A, V, u, b1, Tm, Sm, m, k, tol0=1e-8) {
  for (l in 2:m) { Tm[l, 1:(l-1)][abs(Tm[l, 1:(l-1)]) < tol0] <- 0 }
  for (l in 1:(m-1)) { Tm[l, (l+1):m][abs(Tm[l, (l+1):m]) < tol0] <- 0 }

  # reorder eigenvalues of Schur form!
  eigen_est <- diag(Tm)
  non_zero <- c(which(diag(Tm[-1, ]) != 0))
  for (l in 1:length(non_zero)) {
    eigen_mean <- (eigen_est[non_zero[l]] + eigen_est[non_zero[l] + 1]) / 2
    eigen_est[non_zero[l]] <- eigen_mean
    eigen_est[non_zero[l] + 1] <- eigen_mean
  }
  od <- order(-eigen_est)
  while (od[k] %in% non_zero) {
    eigen_est[od[k]] <- 0
    eigen_est[od[k] + 1] <- 0
    od <- order(-eigen_est)
  }
  od <- order(od)
  orderedT <- ordschur(Sm, Tm, od)
  # Sm %*% Tm %*% t(Sm) - (orderedT$U) %*% orderedT$S %*% t(orderedT$U)
  Tm <- orderedT$S
  Sm <- orderedT$U

  # truncate to a Krylov-Schur decomposition of order k
  V11 <- (V %*% Sm)[,1:k]
  T11 <- Tm[1:k, 1:k]
  u11 <- u / norm2(u)
  b11 <- (b1 %*% Sm)[1:k] * norm2(u)

  # extend to a Krylov decomposition of order m
  Unew <- cbind(V11, u11)
  Snew <- rbind(T11, b11)
  # A %*% Unew[,-(k+1)] - Unew %*% Snew

  for (j in (k+2):(m+1)) {
    v <- A %*% Unew[, j - 1]
    w <- t(Unew) %*% v
    v <- v - Unew %*% w
    v_norm2 <- norm2(v)
    Unew <- cbind(Unew, v / v_norm2)
    Snew <- cbind(rbind(Snew, 0), rbind(w, v_norm2))
  }
  for (l in (1:(m-2))) { Snew[l, max(k+2, l+2):m] <- 0 }
  # A %*% Unew[, -(m+1)] - Unew %*% Snew
  # A %*% Unew[, -(m+1)] - Unew[, -(m+1)] %*% Snew[-(m+1), ] - Unew[, m+1] %*% t(Snew[m+1, ])
  return(list(U=Unew[, -(m+1)], S=Snew[-(m+1), ], u=Unew[, m+1], b=Snew[m+1, ]))
  # A %*% U = U %*% S - u %*% t(b)
}

Arnoldi_Krylov <- function(U, S, u, b, m, k) {
  # Reduce the the Krylov Decomposition S to Hessenberg form by Householder matrix

  # the first column
  if (length(which(S[2:m, 1] != 0)) == 1) {
    V <- diag(m)
    V <- rbind(V[1, ], V[3:(k+1), ], V[2, ], V[(k+2):m,])
    S1 <- t(V) %*% S %*% V
  } else {
    non_zero_indice <- (which(S[2:m, 1] != 0))[1]
    em <- rep(0, m-1)
    em[non_zero_indice] <- 1
    v <- S[2:m, 1] - sign(S[non_zero_indice + 1, 1]) * norm2(S[2:m, 1]) * em
    V <- diag(m)
    V[2:m, 2:m] <- diag(m-1) - 2 * v %*% t(v) / norm2(v)^2
    S1 <- t(V) %*% S %*% V
    S1[1, 3:m] <- 0
    S1[3:m, 1] <- 0
  }

  # other columns
  for (j in 2:(m-2)) {
    v <- S1[(j+1):m, j]
    sigma <- sign(S1[j+1, j]) * norm2(v)
    v <- v + sigma * c(1, rep(0, m-j-1))
    V1 <- diag(m)
    V1[(j+1):m, (j+1):m] <- diag(m-j) - 2 * v %*% t(v) / norm2(v)^2
    S1 <- t(V1) %*% S1 %*% V1
    S1[j, min(j+2, m):m] <- 0
    S1[min(j+2, m):m, j] <- 0
    V <- V %*% V1
  }

  U1 <- U %*% V
  b1 <- t(b) %*% V
  # A %*% U1 - U1 %*% S1 - u %*% b1

  return(list(U=U1, S=S1, u=u, b=b1))
}

check_converge <- function(A, V, H, G, n, k, tol=1e-6) {
  od <- order(-diag(H))
  eigenVal <- diag(H)[od]
  eigenVec=G[, od]

  Ritz <- A %*% V %*% eigenVec - t(eigenVal * t(V %*% eigenVec))
  if (max(apply(Ritz[, 1:k], 2, norm2)) < tol*n) {
    return(list(convergence=TRUE, eigenVal=eigenVal, eigenVec=V %*% eigenVec))
  } else {
    return(list(convergence=FALSE, eigenVal=eigenVal, eigenVec=V %*% eigenVec))
  }
}

eigen_KSLE <- function(A, k, max_iter=10000) {
  # check whether A is a square matrix
  if (nrow(A) != ncol(A)) { stop('Input must be a matrix!') }
  n <- nrow(A)
  m <- min(2 * k, n)
  # randomly initialize b
  b <- runif(n, -1, 1)

  ArnoldiHV <- ArnoldiMethod(A, b, V0=NULL, H0=NULL, r0=NULL, n, m, k, First=TRUE)
  SchurD <- SchurDecom(ArnoldiHV$H, m)
  em <- rep(0, m)
  em[m] <- 1
  expandedUS <- truncate_and_expand(A, ArnoldiHV$V, ArnoldiHV$r, em, SchurD$H, SchurD$G, m, k)
  HessUS <- Arnoldi_Krylov(expandedUS$U, expandedUS$S, expandedUS$u, expandedUS$b, m, k)
  SchurD <- SchurDecom(HessUS$S, m)
  converge <- check_converge(A, HessUS$U, SchurD$H, SchurD$G, n, k)

  iter <- 1
  while (iter < max_iter) {
    if (converge$convergence == TRUE) {
      return(list(eigenVal=converge$eigenVal[1:k], eigenVec=converge$eigenVec[, 1:k]))
    } else {
      iter <- iter + 1
      expandedUS <- truncate_and_expand(A, HessUS$U, HessUS$u, HessUS$b, SchurD$H, SchurD$G, m, k)
      HessUS <- Arnoldi_Krylov(expandedUS$U, expandedUS$S, expandedUS$u, expandedUS$b, m, k)
      SchurD <- SchurDecom(HessUS$S, m)
      converge <- check_converge(A, HessUS$U, SchurD$H, SchurD$G, n, k)
    }
  }
  stop('IRAM not convergence!')
}
