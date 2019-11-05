#' @export
#' @importFrom stats cor cutree var
#' @name hclustvar2
#' @title  Hierarchical clustering of variables from a covariance matrix
#' @description Ascendant hierarchical clustering of a set of variables from a covariance/correlation matrix. 
#' @param x a covariance or correlation matrix.
#' @param init an initial partition (a vector of integers indicating the
#' cluster to which each variable is allocated).
#' @return 
#' \item{height }{a set of p-1 non-decreasing real values: the values of the aggregation criterion. } 
#' \item{clusmat }{a p by p matrix with group
#' memberships where each column k corresponds to the elements of the
#' partition in k clusters.} 
#' \item{merge }{a p-1 by 2 matrix.  Row i of \code{merge} describes the merging of clusters at step i of the clustering.
#' If an element j in the row is negative, then observation -j was merged at
#' this stage.  If j is positive then the merge was with the cluster formed at
#' the (earlier) stage j of the algorithm.  Thus negative entries in \code{merge} indicate agglomerations of singletons, and positive entries
#' indicate agglomerations of non-singletons. }
#' @seealso \code{\link{cutreevar}}, \code{\link{plot.hclustvar}},
#' \code{\link{stability}}
#' @keywords cluster multivariate
#' @examples
#' data(decathlon)
#' x <- cor(decathlon[,1:10])
#' tree <- hclustvar2(x, init=NULL)
#' plot(tree, hang = -1, xlab="", sub="")
#'

hclustvar2 <- function(x, init=NULL) {
  cl <- match.call()
  # test if x is a p times p positive definite symetric matrix
  if (is.null(init)) {
    p <- nrow(x)
    init <- 1:p
    labels <- colnames(x)
  } else {
    p <- max(init)
    labels <- paste("cluster", 1:p, sep = "")
  }
  if (p <= 2) stop("The number of variables must be greater than 2.")
  
  MAXVAL <- 1.0e12
  flag <- rep(1, p)                          # active/dead indicator
  a <- rep(0, p-1)                           # left subnode on clustering
  b <- rep(0, p-1)                           # right subnode on clustering
  ia <- rep(0, p-1)                          # R-compatible version of a
  ib <- rep(0, p-1)                          # R-compatible version of b
  lev <- rep(0, p-1)                         # level or criterion values
  card <- rep(1, p)                          # cardinalities
  order <- rep(0, p)                         # R-compatible order for plotting
  
  diss <- matrix(0,p,p)
  #debut <- Sys.time()
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      A <- which(init==i)
      B <- which(init==j)
      diss[i,j] <-clust_diss2(x,A,B)
      diss[j,i] <- diss[i,j]
    }   
  }
  #Sys.time()-debut
  
  nnsnnsdiss <- getnnsvar(diss, flag)           # call to function getnns
  clusmat <- matrix(0, p, p)                 # cluster memberships
  for (i in 1:p) clusmat[i,p] <- i           # init. trivial partition
  
  for (ncl in (p-1):1) {  
    # check for agglomerable pair
    minobs <- -1;  
    mindis <- MAXVAL;
    for (i in 1:p) {
      if (flag[i] == 1) {
        if (nnsnnsdiss$nndiss[i] < mindis) {
          mindis <- nnsnnsdiss$nndiss[i]
          minobs <- i
        }
      }
    }
    # find agglomerands clus1 and clus2, with former < latter
    if (minobs < nnsnnsdiss$nn[minobs]) {
      clus1 <- minobs
      clus2 <- nnsnnsdiss$nn[minobs] }
    if (minobs > nnsnnsdiss$nn[minobs]) {
      clus2 <- minobs
      clus1 <- nnsnnsdiss$nn[minobs] }
    indicescol<-which(clusmat[,ncl+1]==clus1)
    A <- which(init %in% indicescol)
    indicescol<-which(clusmat[,ncl+1]==clus2)
    B <- which(init %in% indicescol)
    matclus1<-c(A,B)
    # So, agglomeration of pair clus1 < clus2 defines cluster ncl
    #------------------------------------ Block for subnode labels 
    a[ncl] <- clus1                       # aine, or left child node
    b[ncl] <- clus2                       # benjamin, or right child node
    # Now build up ia, ib as version of a, b which is R-compliant
    if (card[clus1] == 1) ia[ncl] <- (-clus1)     # singleton
    if (card[clus2] == 1) ib[ncl] <- (-clus2)     # singleton
    if (card[clus1] > 1) {                # left child is non-singleton
      lastind <- 0
      for (i2 in (p-1):(ncl+1)) {        # Must have p-1 >= ncl+1 here
        if (a[i2] == clus1) lastind <- i2    # Only concerns a[i2]
      }
      ia[ncl] <- p - lastind             # label of non-singleton
    }
    if (card[clus2] > 1) {                # right child is non-singleton
      lastind <- 0
      for (i2 in (p-1):(ncl+1)) {        # Must have p-1 >= ncl+1 here
        if (a[i2] == clus2) lastind <- i2    # Can only concern a[i2]
      }
      ib[ncl] <- p - lastind             # label of non-singleton
    }
    if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
      left <- min(ia[ncl],ib[ncl])
      right <- max(ia[ncl],ib[ncl])
      ia[ncl] <- left                    # Just get left < right
      ib[ncl] <- right
    }
    lev[ncl] <- mindis
    for (i in 1:p) {
      clusmat[i,ncl] <- clusmat[i,ncl+1]
      if (clusmat[i,ncl] == clus2) clusmat[i,ncl] <- clus1
    }
    # Next we need to update diss array
    for (i in 1:p) {
      if ( (i != clus1) && (i != clus2) && (flag[i] == 1) ) {
        indicescol <- which(clusmat[,ncl+1]==i)
        A <- which(init %in% indicescol)
        diss[clus1,i] <- clust_diss2(x, matclus1,A)
        diss[i,clus1] <- diss[clus1,i]
      }
    }
    card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
    # Cluster label clus2 is knocked out; following not nec. but no harm
    flag[clus2] <- 0
    nnsnnsdiss$nndiss[clus2] <- MAXVAL
    for (i in 1:p) {
      diss[clus2,i] <- MAXVAL
      diss[i,clus2] <- diss[clus2,i] }
    # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
    # i.e. nearest neighbors and the nearest neigh. dissimilarity
    nnsnnsdiss <- getnnsvar(diss, flag)
  }
  temp <- cbind(a,b)
  merge2 <- temp[nrow(temp):1, ]
  temp <- cbind(ia,ib)
  merge <- temp[nrow(temp):1,]
  dimnames(merge) <- NULL
  # merge is R-compliant; later suppress merge2
  
  #-------------------------------- Build R-compatible order from ia, ib
  orderlist <- c(merge[p-1,1], merge[p-1,2])
  norderlist <- 2
  for (i in 1:(p-2)) {           # For precisely p-2 further node expansions
    for (i2 in 1:norderlist) {       # Scan orderlist
      if (orderlist[i2] > 0) {     # Non-singleton to be expanded
        tobeexp <- orderlist[i2]
        if (i2 == 1) {
          orderlist <- c(merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[2:norderlist]) }
        if (i2 == norderlist) {
          orderlist <- c(orderlist[1:(norderlist-1)],
                         merge[tobeexp,1],merge[tobeexp,2]) }
        if (i2 > 1 && i2 < norderlist) {
          orderlist <- c(orderlist[1:(i2-1)], 
                         merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[(i2+1):norderlist]) }
        norderlist <- length(orderlist)
      }
    }
  }
  
  orderlist <- (-orderlist)
  class(orderlist) <- "integer"
  xcall <- "hierclust(X)"
  class(xcall) <- "call"
  
  retlist <- list(call = cl, init=init, merge=merge,height=lev[(p-1):1],
                  order=orderlist, labels=labels, clusmat=clusmat)
  class(retlist) <- c("hclust")
  retlist
}
