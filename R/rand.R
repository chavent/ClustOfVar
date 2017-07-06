#' @export
#' @importFrom utils combn
#' @name rand
#' @title Rand index between two partitions
#' @description Returns the Rand index, the corrected Rand index or the asymmetrical Rand
#' index.  The asymmetrical Rand index (corrected or not) measures the
#' inclusion of a partition P into and partition Q with the number of
#' clusters in P greater than the number of clusters in Q.
#' @param P a factor, e.g., the first partition.
#' @param Q a factor, e.g., the second partition.
#' @param symmetric a boolean. If FALSE the asymmetrical Rand index is
#' calculated.
#' @param adj a boolean. If TRUE the corrected index is calculated.
#' @seealso \code{\link{stability}}

rand <- function(P,Q,symmetric=TRUE,adj=TRUE)
{
  tab.cont <- table(P,Q)
  cptnij <- 0
  for (i in 1:nrow(tab.cont)) {
    for   (j in 1:ncol(tab.cont))  {
      cptnij <- cptnij+(tab.cont[i,j]*(tab.cont[i,j]-1))/2
    }
  }
  cptni <- 0
  cptnj <- 0
  ni. <- apply(tab.cont,1,sum)
  n.j <- apply(tab.cont,2,sum)
  for (i in 1:nrow(tab.cont)) {
    cptni <- cptni+(ni.[[i]]*(ni.[[i]]-1))/2
  }
  for (j in 1:ncol(tab.cont)) {
    cptnj <- cptnj+(n.j[[j]]*(n.j[[j]]-1))/2
  }
  a <- cptnij
  b <- cptni-cptnij
  cc <- cptnj-cptnij
  n <- sum(ni.)
  d <- (n*(n-1))/2+cptnij-cptni-cptnj
  
  if (symmetric==TRUE & adj==FALSE){
    return((a+d)/(a+b+cc+d))
  }
  
  if (symmetric==FALSE & adj==FALSE){
    return((a+d+cc)/(a+b+cc+d))
  }
  
  if (symmetric==TRUE & adj==TRUE){
    return((cptnij-((cptni*cptnj)/ncol(combn(n,2))))/(1/2*(cptni+cptnj)-((cptni*cptnj)/((n*(n-1))/2))))
  }
  
  if (symmetric==FALSE & adj==TRUE){
    return((cptnij-((cptni*cptnj)/ncol(combn(n,2))))/(cptni-((cptni*cptnj)/((n*(n-1))/2))))
  }
}

