##' @export
##' @method plot hclustvar
##' @title Dendrogram of the hierarchy of variables 
##' @description Dendrogram of the hierarchy of variables resulting from \code{hclustvar}
##' and aggregation levels plot.
##' 
##' 
##' @param x an object of class \code{hclustvar}.
##' @param type if type="tree" plot of the dendrogram and if type="index"
##' aggregation levels plot.
##' @param which if one of the two plots is required, specify a subset of the
##' numbers 1:2
##' @param ask logical; if TRUE, the user is _ask_ed before each plot.
##' @param sub a sub title for the plot.
##' @param \dots further arguments passed to or from other methods.
##' @seealso \code{\link{hclustvar}}
##' @keywords cluster multivariate
##' @examples
##' 
##' data(wine)
##' tree <- hclustvar(data=wine)
##' plot(tree)
##' 
##' # 2 plots on 1 page
##' par(mfrow = c(1, 2))
##' plot(tree) 
##' 
##'  # plot just the dendrogram
##' plot(tree,which=1)
##' 
plot.hclustvar <-
  function(x,type="tree",which=c(1:2), ask = prod(par("mfcol")) < length(which) && dev.interactive(),sub="", ...){
    if (!inherits(x, "hclustvar")) 
      stop("use only with \"hclustvar\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 2)) 
      stop("'which' must be in 1:2")
    class(x) <- "hclust"
    #show <- rep(FALSE, 2)
    #show[which] <- TRUE
    #if (ask) {
    #  oask <- devAskNewPage(TRUE)
    #  on.exit(devAskNewPage(oask))
    # }	
    #if (show[1]) plot(x, hang=-1, xlab="", sub=sub, ...)	
    #if (show[2]) {
    #	plot(x=seq(length(x$height),1),x$height,xaxt = "n",ylim=c(0,max(x$height)),xlab="number of clusters",ylab="Height",main="Aggregation levels",type="n")
    #	points(x=seq(length(x$height),1),x$height,pch=3)
    #	axis(side=1,at=seq(1,length(x$height)),labels=paste(1:(length(x$height))))
    #}
    if ((type !="tree") && (type !="index")) { stop("type must be equal to  \"tree \" or  \"index\"")} else {
      if (type=="tree") plot(x, hang=-1, xlab="", sub=sub, ...) 
      if (type=="index"){
        plot(x=seq(length(x$height),1),x$height,xaxt = "n",ylim=c(0,max(x$height)),xlab="number of clusters",ylab="Height",main="Aggregation levels",type="n")
        points(x=seq(length(x$height),1),x$height,pch=3)
        axis(side=1,at=seq(1,length(x$height)),labels=paste(1:(length(x$height)))) 
       }
    }
    
    class(x) <- c("hclustvar","hclust")
  }
