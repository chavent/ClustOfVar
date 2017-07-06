#' @export
#' @name stability
#' @title Stability of partitions from a hierarchy of variables
#' @description Evaluates the stability of partitions obtained from a hierarchy of \code{p}
#' variables. This hierarchy is performed with \code{hclustvar} and the
#' stability of the partitions of 2 to p-1 clusters is evaluated with a
#' bootstrap approach. The boostrap approch is the following: \code{hclustvar}
#' is applied to \code{B} boostrap samples of the \code{n} rows.  The
#' partitions of 2 to p-1 clusters obtained from the B bootstrap hierarchies
#' are compared with the partitions from the initial hierarchy . The mean of
#' the corrected Rand indices is plotted according to the number of clusters.
#' This graphical representation helps in the determination of a suitable
#' numbers of clusters.
#' @param tree an object of class \code{hclustvar}.
#' @param B the number of bootstrap samples.
#' @param graph boolean, if 'TRUE' a graph is displayed.
#' @return 
#' \item{matCR}{matrix of corrected Rand indices.}
#' \item{meanCR}{vector of mean corrected Rand indices.}
#' @seealso \code{\link{plot.clustab}}, \code{\link{hclustvar}}
#' @examples
#' data(decathlon)
#' tree <- hclustvar(X.quanti=decathlon[,1:10])
#' stab<-stability(tree,B=20)
#' plot(stab,nmax=7)
#' boxplot(stab$matCR[,1:7])
#'
stability <- function(tree,B=100,graph=TRUE)
  {
    if (!inherits(tree, "hclustvar")) 
      stop("use only with \"hclustvar\" objects")
    cl <- match.call()
    X.quali <- tree$X.quali
    X.quanti <- tree$X.quanti
    clusmat <- tree$clusmat
    nmax <- ncol(clusmat)-1
    matRandCorrige <- matrix(NA,nrow=B,ncol=nmax-1)
    matRandNonCorrige <- matrix(NA,nrow=B,ncol=nmax-1)
    for (b in 1:B)
    {
      #print(b)
      res1<-bootvar(X.quanti,X.quali)
      Xboot.quanti<-res1$Xboot.quanti
      Xboot.quali<-res1$Xboot.quali
      # a ameliorer peut-etre
      if (!is.null(Xboot.quanti)&& !is.null(Xboot.quali)) 
      {clusmatboot <- hclustvar(Xboot.quanti,Xboot.quali)$clusmat}
      if (is.null(Xboot.quanti)&& !is.null(Xboot.quali)) 
      {clusmatboot <- hclustvar(X.quali =Xboot.quali)$clusmat}
      if (!is.null(Xboot.quanti)&& is.null(Xboot.quali)) 
      {clusmatboot <- hclustvar(X.quanti=Xboot.quanti)$clusmat}	
      for (i in 2:nmax) 
      {matRandCorrige[b,i-1] <- rand(clusmat[,i],clusmatboot[,i],adj=TRUE)}
      for (i in 2:nmax) 
      {matRandNonCorrige[b,i-1] <- rand(clusmat[,i],clusmatboot[,i],adj=FALSE)}
    }
    colnames <- paste("P", 2:nmax,sep = "")
    rownames <- paste("B", 1:B, sep = "")
    colnames(matRandCorrige) <- colnames
    rownames(matRandCorrige) <- rownames
    critRandmoyenCorrige<-apply(matRandCorrige,2,mean)
    if (graph!=FALSE) {
      plot(x=seq(1,nmax-1),critRandmoyenCorrige,xaxt="n",ylim=c(0,1),xlab="number of clusters",ylab="mean adjusted Rand criterion",main="Stability of the partitions",type="b")
      axis(side=1,at=seq(1,nmax-1),labels=paste(2:nmax))
    }

    retlist <- list(call=cl,matCR=matRandCorrige,meanCR=critRandmoyenCorrige)
    class(retlist) <- "clustab"
    retlist
  }

