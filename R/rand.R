##' @export
##' @name rand
##' @title Rand index between two partitions
##' 
##' @description Returns the Rand index, the corrected Rand index or the asymmetrical Rand
##' index.  The asymmetrical Rand index (corrected or not) measures the
##' inclusion of a partition 'P' into and partition 'Q' with the number of
##' clusters in 'P' greater than the number of clusters in 'Q'.
##' 
##' 
##' @param P a factor, e.g., the first partition.
##' @param Q a factor, e.g., the second partition.
##' @param symmetric a boolean. If 'FALSE' the asymmetrical Rand index is
##' calculated.
##' @param adj a boolean. If 'TRUE' the corrected index is calculated.
##' @author Marie Chavent <marie.chavent@@u-bordeaux.fr>, Vanessa Kuentz, Amaury Labenne,
##' Benoit Liquet, Jerome Saracco
##' @seealso \code{\link{stability}}
##' @keywords cluster multivariate
rand <-
function(P,Q,symmetric=TRUE,adj=TRUE)
{
	#calculs communs aux versions symétriques et non symétriques
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
	#Rand symetrique non corrigé
	if (symmetric==TRUE & adj==FALSE){
  		return((a+d)/(a+b+cc+d))
	}
	#Rand asymetrique non corrigé
	if (symmetric==FALSE & adj==FALSE){
  		return((a+d+cc)/(a+b+cc+d))
	}
	#Rand symetrique corrigé
	if (symmetric==TRUE & adj==TRUE){
  		return((cptnij-((cptni*cptnj)/ncol(combn(n,2))))/(1/2*(cptni+cptnj)-((cptni*cptnj)/((n*(n-1))/2))))
	}
	#Rand asymetrique corrigé
	if (symmetric==FALSE & adj==TRUE){
  		return((cptnij-((cptni*cptnj)/ncol(combn(n,2))))/(cptni-((cptni*cptnj)/((n*(n-1))/2))))
	}
}

