#' Call Schur() to decompose and sort
#' @param H input matrix
#' @return A list containing Q and T for H = QTQ'
#' @examples
#' H <- matrix(c(1,2,3,4,2,5,6,7,3,6,8,9,4,7,9,10), 4, 4)
#' SchurDecom(H)
#' @import Matrix
#' @export
#'
SchurDecom <- function(H) {
  SchurOut <- Schur(H)
  Q1 <- SchurOut$Q
  T1 <- SchurOut$T

  # reorder
  od <- order(-diag(T1))
  T1 <- diag(diag(T1)[od])
  Q1 <- Q1[, od]
  return(list(Q=Q1, T=T1))
}
