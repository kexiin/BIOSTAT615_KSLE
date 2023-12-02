truncatedSVD <- function(M, k, method='KSLE') {
  if (method == 'IRAM') {
    UD <- eigen_IRAM(M %*% t(M), k)
    VD <- eigen_IRAM(t(M) %*% M, k)
    return(list(D=sqrt(UD$eigenVal), U=UD$eigenVec, V=VD$eigenVec))
  } else {
    UD <- eigen_KSLE(M %*% t(M), k)
    VD <- eigen_KSLE(t(M) %*% M, k)
    return(list(D=sqrt(UD$eigenVal), U=UD$eigenVec, V=VD$eigenVec))
  }
}
