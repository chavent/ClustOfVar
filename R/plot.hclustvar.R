#' @export
#' @name plot.hclustvar
#' @importFrom graphics axis dotchart plot points
#' @title Dendrogram of the hierarchy of variables
#' @description Dendrogram of the hierarchy of variables resulting from \code{hclustvar}
#' and aggregation levels plot.
#' @param x an object of class \code{hclustvar}.
#' @param type if type="tree" plot of the dendrogram and if type="index"
#' aggregation levels plot.
#' @param sub a sub title for the plot.
#' @param \dots further arguments passed to or from other methods.
#' @seealso \code{\link{hclustvar}}
#' @examples
#' data(wine)
#' X.quanti <- PCAmixdata::splitmix(wine)$X.quanti
#' X.quali <- PCAmixdata::splitmix(wine)$X.quali
#' tree <- hclustvar(X.quanti,X.quali)
#' plot(tree)
#'
#' #Aggregation levels plot
#' plot(tree,type="index")
#'
plot.hclustvar <- function(x,type="tree",sub="", ...){
  if (!inherits(x, "hclustvar")) 
    stop("use only with \"hclustvar\" objects")
  class(x) <- "hclust"
  
  if (type=="tree") plot(x, hang=-1, xlab="", sub=sub, ...) 
  if (type=="index"){
    plot(x=seq(length(x$height),1),x$height,xaxt = "n",ylim=c(0,max(x$height)),xlab="number of clusters",ylab="Height",main="Aggregation levels",type="n")
    points(x=seq(length(x$height),1),x$height,pch=3)
    axis(side=1,at=seq(1,length(x$height)),labels=paste(1:(length(x$height)))) 
  }
  
  class(x) <- c("hclustvar","hclust")
}
