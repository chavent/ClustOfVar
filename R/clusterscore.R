#' @export
#' @name clusterscore
#' @title Calculates de synthetic variable of a cluster 
#' @description Calculates the synthetic variable of a cluster of variables.
#' The variables can be quantitative or qualitative. The synthetic variable is the first principal component of PCAmix.
#' The variance of the synthetic variable is the first eigenvalue. It is equal to the sum of 
#' squared correlations or correlation ratios to the synthetic variable. It measures the homogeneity of the cluster.
#' @param Z a centered and reduced data matrix obtained with the recod function 
#' @return \item{f}{the synthetic variables i.e. the scores on the first principal component of PCAmix}
#' @return \item{sv}{the standard deviation of f i.e. the first singular value}
#' @return \item{v}{the standardized loadings}
#' @examples 
#' data(decathlon)
#' A <- 1:5
#' Z <- PCAmixdata::recod(X.quanti=decathlon[1:10,A], X.quali=NULL)$Z
#' clusterscore(Z)
#' Z%*%as.matrix(clusterscore(Z)$v)
#' clusterscore(Z)$f

clusterscore <- function(Z)
{
  
  Z<-as.matrix(Z)
  n<-nrow(Z)
  Ztilde <-Z/sqrt(n)
  e <- svd(Ztilde)
  f <- e$u[,1]*e$d[1]*sqrt(n) 
  sv <-e$d[1]
  v <- e$v[,1]
  return(list(f=f,sv=sv,v=v))	
}

