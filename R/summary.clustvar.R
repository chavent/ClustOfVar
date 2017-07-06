#' @export summary.clustvar
#' @export
#' @name summary.clustvar
#' @title Summary of a 'hclustvar' object
#' @description This is a method for the function summary for objects of the class
#' \code{clustvar}.
#' @param object an object of class \code{clustvar}.
#' @param \dots further arguments passed to or from other methods.
#' @return Returns a list of matrices of squared loadings i.e. for each
#' cluster of variables, the squared loadings on first principal component of
#' PCAmix. For quantitative variables (resp. qualitative), squared loadings
#' are the squared correlations (resp. the correlation ratios) with the first
#' PC (the cluster center). If the partition of variables has been obtained
#' with \code{kmeansvar} the number of iteration until convergence is also
#' indicated.
#' @seealso \code{\link{kmeansvar}}, \code{\link{cutreevar}}
#' @examples
#' data(decathlon)
#' part<-kmeansvar(X.quanti=decathlon[,1:10],init=5)
#' summary(part)
#'
summary.clustvar <- function(object, ...)
{
  x <- object
  if (!inherits(object, "clustvar")) 
    stop("use only with \"clustvar\" objects")
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  iter <- object$iter
  if (is.numeric(iter)) cat("number of iterations: ",iter,sep=" ") 
  cat("\n")
  cat("\n")
  n <- object$rec$n
  p <- object$rec$p
  p1 <- object$rec$p1
  p2 <- p-p1
  cat("Data:", "\n")
  cat(paste("   number of observations: ",n),sep=" ") 
  cat("\n")
  if 	((p1!=0)&& (p2==0)) {
    cat(paste("   number of variables: ",p1),sep=" ")  
    cat("\n")
  }
  if 	((p1==0)&& (p2!=0)) {
    cat(paste("   number of variables: ",p2),sep=" ")  
    cat("\n")
  }
  if 	((p1!=0)&& (p2!=0)) {
    cat(paste("   number of  variables: ",p),sep=" ")
    cat("\n")
    cat(paste("        number of numerical variables: ",p1),sep=" ")   
    cat("\n")
    cat(paste("        number of categorical variables: ",p2),sep=" ")   
    cat("\n")
  }
  cat(paste("   number of clusters: ",object$k),sep=" ") 
  cat("\n")
  
  
  for (g in 1:object$k)
  {
    cat("\n")
    cat(paste("Cluster ",g,": "),sep=" ")
    cat("\n")
    print(object$var[[g]],digits=2)
    cat("\n")
  }
  cat("\n")
  cat(paste("Gain in cohesion (in %): ",round(object$E,digits=2)),sep=" ") 
  cat("\n")
}

