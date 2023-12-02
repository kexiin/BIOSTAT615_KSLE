norm2 <- function(x) {
  return(drop(sqrt(sum(x^2))))
}

ArnoldiMethod <- function(A, b, V0, H0, r0, n, m, k, First=FALSE) {
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
  for (j in start:(m+1)) {
    Vj <- A %*% V[, j-1]
    for (i in max(1, j-2):(j-1)) {
      H[i, j-1] <- t(V[, i]) %*% Vj
      Vj <- Vj - H[i, j-1] * V[, i]
    }
    H[j, j-1] <- norm2(Vj)
    V[, j] <- Vj / H[j, j-1]
  }
  r <- V[, m+1] * H[m+1, m]
  V <- V[, -(m+1)]
  H <- H[-(m+1), ]
  for (i in 3:m) { H[i, 1:(i-2)] <- 0 }
  for (i in 1:(m-2)) { H[i, (i+2):m] <- 0 }
  # A %*% V - V %*% H - r %*% t(e)

  return(list(H=H, V=V, r=r))
}

Arnoldi_converge <- function(A, b, V0, H0, r0, n, m, k, First=FALSE, tol=1e-6) {
  ArnoldiVH <- ArnoldiMethod(A, b, V0, H0, r0, n, m, k, First)
  H <- ArnoldiVH$H
  V <- ArnoldiVH$V
  r <- ArnoldiVH$r

  # Solve the eigenpairs of Hessenberg Matrix H
  eigenpairsH <- shiftedQR(H, m)
  # Ritz Vector conditions for termination
  Ritz <- A %*% V %*% eigenpairsH$eigenVec - t(eigenpairsH$eigenVal * t(V %*% eigenpairsH$eigenVec))
  if (max(apply(Ritz[, 1:k], 2, norm2)) < tol*n) {
    return(list(convergence=TRUE, eigenVal=eigenpairsH$eigenVal, eigenVec=V %*% eigenpairsH$eigenVec))
  } else {
    return(list(convergence=FALSE, eigenVal=eigenpairsH$eigenVal, H=H, V=V, r=r, ritz=Ritz[,1:k]))
  }
}


GivensRotation <- function(H1, sigma, m) {
  H2 <- H1 - sigma * diag(m)
  G <- diag(m)
  for (i in 1:(m-1)) {
    r <- sqrt((H2[i, i])^2 + (H2[i + 1, i])^2)
    c <- H2[i, i] / r
    s <- H2[i + 1, i] / r
    Gi <- diag(m)
    Gi[i:(i+1), i:(i+1)] <- matrix(c(c, -s, s, c), 2, 2)
    H2 <- Gi %*% H2
    G <- Gi %*% G
  }
  return(list(H_new=H2 %*% t(G) + sigma * diag(m), G=G))
}

shiftedQR <- function(H, m, tol=1e-6, tol0=1e-14, max_iter=10000, check='weak') {
  # initialize
  Hlast <- diag(m)
  H1 <- H
  iter <- 0
  G1 <- diag(m)
  changeH1 <- ifelse(check=='weak', ((norm2(diag(H1) - diag(Hlast))) > tol),
                     ((norm2(H1 - Hlast) > tol*1e3) | ((norm2(diag(H1) - diag(Hlast))) > tol)))
  # changeH1 <- ((norm2(diag(H1) - diag(Hlast))) > tol)

  # H = G1 %*% H1 %*% t(G1), where H1 is upper tri, G1'G1=I
  while (changeH1) {
    iter <- iter + 1
    Hlast <- H1
    if (iter > max_iter) {
      stop('Shifted QR not convergence!')
    }

    # Wilkinson shift
    non_zero <- c(which(abs(diag(H1[-1, ])) > tol0))
    shift_m <- ifelse(length(non_zero) > 0, non_zero[length(non_zero)] + 1, m)
    eta <- (H1[shift_m-1, shift_m-1] - H1[shift_m, shift_m]) / 2
    sigma <- H1[shift_m, shift_m] + eta - sign(eta) * sqrt(eta^2 + (H1[shift_m-1,shift_m])^2)
    # sigma <- H1[m, m]

    # Givens rotation
    qrH <- GivensRotation(H1, sigma, m)
    H1 <- qrH$H_new
    G1 <- G1 %*% t(qrH$G)
    changeH1 <- ifelse(check=='weak', ((norm2(diag(H1) - diag(Hlast))) > tol),
                       ((norm2(H1 - Hlast) > tol*1e3) | ((norm2(diag(H1) - diag(Hlast))) > tol)))
    # changeH1 <- ((norm2(diag(H1) - diag(Hlast))) > tol)
  }

  od <- order(-diag(H1))
  eigenVal <- diag(H1)[od]
  return(list(eigenVal=eigenVal, eigenVec=G1[, od], H=H1, G=G1))
}

ImplicitlyRestart <- function(H, V, r, eigenVal, k, m) {
  H1 <- H
  V1 <- V
  G <- diag(m)

  # restart by removing unwanted eigenvalues
  for (i in (k+1):m) {
    sigma <- eigenVal[i]
    qrH <- GivensRotation(H1, sigma, m)
    H1 <- qrH$H_new
    for (i in 3:m) { H1[i, 1:(i-2)] <- 0 }
    for (i in 1:(m-2)) { H1[i, (i+2):m] <- 0 }
    V1 <- V1 %*% t(qrH$G)
    G <- G %*% t(qrH$G)
  }
  r1 <- V1[, k+1] * H1[k+1, k] + r * G[m, k]
  # A %*% V1[, 1:k] - V1[,1:k] %*% H1[1:k, 1:k] - r1 %*% t(e)

  return(list(H=H1[1:k, 1:k], V=V1[, 1:k], r=r1))
}


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
      return(list(eigenVal=eigenH$eigenVal[1:k], eigenVec=eigenH$eigenVec[, 1:k]))
    } else {
      iter <- iter + 1
      newVHr <- ImplicitlyRestart(eigenH$H, eigenH$V, eigenH$r, eigenH$eigenVal, k, m)
      eigenH <- Arnoldi_converge(A, b, newVHr$V, newVHr$H, newVHr$r, n, m, k, First=FALSE)
    }
  }
  stop('IRAM not convergence!')
}
