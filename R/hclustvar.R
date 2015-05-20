##' @export
##' @name hclustvar
##' @title  Hierarchical clustering of variables
##' 
##' @description Ascendant hierarchical clustering of a set of variables.  Variables can be
##' quantitative, qualitative or a mixture of both. The aggregation criterion
##' is the decrease in homogeneity for the clusters being merged. The
##' homogeneity of a cluster is the sum of the correlation ratio (for
##' qualitative variables) and the squared correlation (for quantitative
##' variables) between the variables and the center of the cluster which is the
##' first principal component of PCAmix. PCAmix is defined for a mixture of
##' qualitative and quantitative variables and includes ordinary principal
##' component analysis (PCA) and multiple correspondence analysis (MCA) as
##' special cases. Missing values are replaced by means for quantitative
##' variables and by zeros in the indicator matrix for qualitative variables.
##' 
##' 
##' @param X.quanti a numeric matrix of data, or an object that can be coerced
##' to such a matrix (such as a numeric vector or a data frame with all numeric
##' columns).
##' @param X.quali a categorical matrix of data, or an object that can be
##' coerced to such a matrix (such as a character vector, a factor or a data
##' frame with all factor columns).
##' @param init an initial partition (a vector of integers indicating the
##' cluster to which each variable is allocated).
##' @return \item{height }{a set of p-1 non-decreasing real values: the values
##' of the aggregation criterion. } \item{clusmat }{a p by p matrix with group
##' memberships where each column k corresponds to the elements of the
##' partition in k clusters.} \item{merge }{a p-1 by 2 matrix.  Row i of
##' \code{merge} describes the merging of clusters at step i of the clustering.
##' If an element j in the row is negative, then observation -j was merged at
##' this stage.  If j is positive then the merge was with the cluster formed at
##' the (earlier) stage j of the algorithm.  Thus negative entries in
##' \code{merge} indicate agglomerations of singletons, and positive entries
##' indicate agglomerations of non-singletons. }
##' @details If the quantitative and qualitative data are in a same dataframe, the function 
##' \code{splitmix} can be used to extract automatically the qualitative and the quantitative 
##' data in two separated dataframes.
##' @author Marie Chavent \email{marie.chavent@@u-bordeaux.fr}, Vanessa Kuentz, Amaury Labenne
##' Benoit Liquet, Jerome Saracco
##' @seealso \code{\link{cutreevar}}, \code{\link{plot.hclustvar}},
##' \code{\link{stability}},\code{\link{kmeansvar}},\code{\link{splitmix}}
##' @references Chavent, M., Liquet, B., Kuentz, V., Saracco, J. (2012),
##' ClustOfVar: An R Package for the Clustering of Variables. Journal of
##' Statistical Software, Vol. 50, pp. 1-16.
##' 
##' Chavent, M., Kuentz, V., Saracco, J. (2011), Orthogonal Rotation in PCAMIX.
##' Advances in Classification and Data Analysis, Vol. 6, pp. 131-146.
##' @keywords cluster multivariate
##' @examples 
##' 
##' #quantitative variables
##' data(decathlon)
##' tree <- hclustvar(X.quanti=decathlon[,1:10], X.quali=NULL, init=NULL)
##' plot(tree)
##' 
##' #qualitative variables with missing values
##' data(vnf) 
##' tree_NA <- hclustvar(X.quali=vnf, X.quanti=NULL)
##' plot(tree_NA)
##' dev.new()
##' vnf2<-na.omit(vnf)
##' tree <- hclustvar(X.quali=vnf2, X.quanti=NULL)
##' plot(tree)
##' 
##' #mixture of quantitative and qualitative variables
##' data(wine)
##' X.quanti <- splitmix(wine)$X.quanti[,1:27] 
##' X.quali <- splitmix(wine)$X.quali
##' tree <- hclustvar(X.quanti,X.quali)
##' plot(tree)
##' 
##' #combined clustering
##' #library(mixOmics)
##' #data(breast.tumors)
##' #X.quanti <- breast.tumors$gene.exp
##' #init<- kmeansvar(X.quanti,init=100)$cluster
##' #tree <- hclustvar(X.quanti,init=init)
##' #plot(tree)
##' 
##' #data(yeast)
##' #X.quanti <- yeast$data
##' #tree <- hclustvar(X.quanti)
##' #plot(tree)
##' #init<- cutreevar(tree,k=10)$cluster
##' #tree2 <- hclustvar(X.quanti,init=init)
##' #plot(tree2)
##' 
##' #cutreevar(tree,3)$cluster
##' #cutreevar(tree2,3)$cluster
##' @import PCAmixdata
hclustvar <- function(X.quanti=NULL,X.quali=NULL,init=NULL) {
  cl <- match.call()
  rec <- recod(X.quanti,X.quali,rename.level=TRUE)
  X <- rec$X     
  Z <- rec$Z		
  indexj <- rec$indexj 
  n <- rec$n
  if (is.null(init)) {
    p <- rec$p
    init <- 1:p
    labels <- colnames(X)
  } else {
    p <- max(init)
    labels <- paste("cluster", 1:p, sep = "")
  }
 
  #X.quanti<-rec$X.quanti
  #X.quali<-rec$X.quali
  
  
  if (p<=2) stop("The number of variables must be greater than 2.")
  MAXVAL <- 1.0e12
  flag <- rep(1, p)                          # active/dead indicator
  a <- rep(0, p-1)                           # left subnode on clustering
  b <- rep(0, p-1)                           # right subnode on clustering
  ia <- rep(0, p-1)                          # R-compatible version of a
  ib <- rep(0, p-1)                          # R-compatible version of b
  lev <- rep(0, p-1)                         # level or criterion values
  card <- rep(1, p)                          # cardinalities
  order <- rep(0, p)                         # R-compatible order for plotting
  
  diss<-matrix(0,p,p)
  debut <- Sys.time()
  for (i in 1:(p-1)) {
    a <- i+1
    for (j in a:p) {
      A <- which(init==i)
      B <- which(init==j)
      matA<-Z[,which(is.element(indexj,A))]
      matB <- Z[,which(is.element(indexj,B))]
      diss[i,j] <-clust_diss(matA,matB)
      diss[j,i] <- diss[i,j]
    }   
  }
  Sys.time()-debut
  
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
    Xclus1<-Z[,which(is.element(indexj,A))]
    indicescol<-which(clusmat[,ncl+1]==clus2)
    B <- which(init %in% indicescol)
    Xclus2<-Z[,which(is.element(indexj,B))]
    matclus1<-cbind(Xclus1,Xclus2)
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
        mati<-Z[,which(is.element(indexj,A))]
        diss[clus1,i] <- clust_diss(matclus1,mati)
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
  
  retlist <- list(call = cl,rec=rec,init=init,merge=merge,height=lev[(p-1):1],order=orderlist,labels=labels,clusmat=clusmat,X.quanti=X.quanti,X.quali=X.quali)
  class(retlist) <- c("hclustvar","hclust")
  retlist
}
