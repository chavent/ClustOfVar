#' @export
#' @name clust_diss
#' @title Calculates the aggregation criterion for two clusters of variables
#' @description Calculates the measure of aggregation of two clusters of variables. This
#' measure of aggregation is equal to the decrease in homogeneity for the
#' clusters being merged.
#' @param A a centered and reduced data matrix obtained with the recod
#' function for the first cluster
#' @param B a centered and reduced data matrix obtained with the recod
#' function for the second cluster
#' @return The aggregation measure between the two clusters
#' @examples
#' data(decathlon)
#' A <- PCAmixdata::recod(X.quanti=decathlon[1:10,1:5], X.quali=NULL)$Z
#' B <- PCAmixdata::recod(X.quanti=decathlon[1:10,6:10], X.quali=NULL)$Z
#' clust_diss(A,B)

clust_diss <- function(A,B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  AUB <- cbind(A,B)
  repA <- clusterscore(A)
  repB <- clusterscore(B)
  repAUB <- clusterscore(AUB)
  valproA <- repA$sv^2
  valproB <- repB$sv^2
  valproAUB <- repAUB$sv^2
  crit <- valproA+valproB-valproAUB
  if (crit < 1e-7) crit <- 0
  return(crit)
}

