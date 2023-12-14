#' A implementation of SVD, both IRAM and Krylov-Schur algorithm
#' @param M input matrix
#' @param k need k singular values
#' @param method choice includes "IRAM" and "KSLE"
#' @param Uonly whether only need U and D
#' @return A list containing U, D and V
#' @export
#'
truncatedSVD <- function(M, k, method='KSLE', Uonly=FALSE) {
  if (Uonly == FALSE) {
    if (method == 'IRAM') {
      UD <- eigen_IRAM(M %*% t(M), k)
      VD <- eigen_IRAM(t(M) %*% M, k)
      return(list(D=sqrt(UD$eigenVal), U=UD$eigenVec, V=VD$eigenVec))
    } else {
      UD <- eigen_KSLE(M %*% t(M), k)
      VD <- eigen_KSLE(t(M) %*% M, k)
      return(list(D=sqrt(UD$eigenVal), U=UD$eigenVec, V=VD$eigenVec))
    }
  } else {
    if (method == 'IRAM') {
      UD <- eigen_IRAM(M %*% t(M), k)
      return(list(D=sqrt(UD$eigenVal), U=UD$eigenVec))
    } else {
      UD <- eigen_KSLE(M %*% t(M), k)
      return(list(D=sqrt(UD$eigenVal), U=UD$eigenVec))
    }
  }
}
