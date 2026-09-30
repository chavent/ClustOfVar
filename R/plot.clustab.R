#' @export
#' @name plot.clustab
#' @title  Plot of an index of stability of partitions of variables
#' @description Plot of the index of stability of the partitions against the number of
#' clusters.
#' @param x an object of class \code{clustab}.
#' @param nmin the minimum number of clusters in the plot.
#' @param nmax the maximum number of clusters in the plot.
#' @param \dots further arguments passed to or from other methods.
#' @seealso \code{\link{stability}}
#' @examples
#' data(decathlon)
#' tree <- hclustvar(X.quanti=decathlon[,1:10])
#' stab<-stability(tree,B=20)
#' plot(stab,nmax=7)
#'
plot.clustab <- function(x,nmin=NULL,nmax=NULL, ...)
{
  if (!inherits(x, "clustab")) 
    stop("use only with \"clustab\" objects")
  if (is.null(nmin)) nmin <- 2
  if (is.null(nmax)) nmax <- length(x$meanCR)+1
  meanCR <- x$meanCR[(nmin-1):(nmax-1)]
  plot(x=seq(nmin,nmax),meanCR,xaxt="n",ylim=c(0,1),xlab="number of clusters",ylab="mean adjusted Rand criterion",type="b",...)
  axis(side=1,at=seq(nmin,nmax),labels=paste(nmin:nmax))
}

