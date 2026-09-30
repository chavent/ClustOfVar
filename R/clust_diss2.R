#' @export
#' @name clust_diss2
#' @title Dissimilarity between two clusters of variables 
#' @description Dissimilarity between two clusters of variables
#'  when only the covariance/correlation matrix is known.
#' @param x a covariance/correlation matrix
#' @param A indices of cluster A
#' @param B indices of cluter B
#' @return The dissimilarity between the two clusters
#' @examples
#' data(decathlon)
#' x <- cor(decathlon[,1:10])
#' A <- c(1,3,4)
#' B <- c(2,7,10)
#' clust_diss2(x,A,B)

clust_diss2 <- function(x, A, B) {
  vpA <- eigen(x[A,A])$values[1]
  vpB <- eigen(x[B,B])$values[1]
  AUB <- c(A,B)
  vpAUB <- eigen(x[AUB,AUB])$values[1]
  crit <- vpA + vpB - vpAUB
  if (crit < 1e-7) crit <- 0
  return(crit)
}

