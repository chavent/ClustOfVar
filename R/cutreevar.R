#' @export
#' @name cutreevar
#' @title Cut a tree into groups of variables
#' @description Cuts a hierarchical tree of variables resulting from \code{hclustvar} into
#' several clusters by specifying the desired number of clusters.
#' @param obj an object of class 'hclustvar'.
#' @param k an integer scalar with the desired number of clusters.
#' @param matsim boolean, if TRUE, the matrices of similarities between
#' variables in same cluster are calculated.
#' @return 
#' \item{var}{a list of matrices of squared loadings i.e. for each
#' cluster of variables, the squared loadings on first principal component of
#' PCAmix. For quantitative variables (resp. qualitative), squared loadings
#' are the squared correlations (resp. the correlation ratios) with the first
#' PC (the cluster center).}
#' \item{sim}{a list of matrices of similarities
#' i.e. for each cluster, similarities between their variables.  The
#' similarity between two variables is defined as a square cosine: the square
#' of the Pearson correlation when the two variables are quantitative; the
#' correlation ratio when one variable is quantitative and the other one is
#' qualitative; the square of the canonical correlation between two sets of
#' dummy variables, when the two variables are qualitative. \code{sim} is
#' 'NULL' if 'matsim' is 'FALSE'. }
#' \item{cluster}{a vector of integers
#' indicating the cluster to which each variable is allocated.} 
#' \item{wss}{the
#' within-cluster sum of squares for each cluster: the sum of the correlation
#' ratio (for qualitative variables) and the squared correlation (for
#' quantitative variables) between the variables and the center of the
#' cluster.} 
#' \item{E}{the pourcentage of homogeneity which is accounted by the
#' partition in k clusters.} 
#' \item{size}{the number of variables in each
#' cluster.} 
#' \item{scores}{a n by k numerical matrix which contains the k
#' cluster centers. The center of a cluster is a synthetic variable: the first
#' principal component calculated by PCAmix.  The k columns of \code{scores}
#' contain the scores of the n observations units on the first PCs of the k
#' clusters.} 
#' \item{coef}{a list of the coefficients of the linear
#' combinations defining the synthetic variable of each cluster.}
#' @seealso \code{\link{hclustvar}}
#' @examples
#' data(decathlon)
#' tree <- hclustvar(decathlon[,1:10])
#' plot(tree)
#' #choice of the number of clusters
#' stability(tree,B=40)
#' part <- cutreevar(tree,4)
#' print(part)
#' summary(part)
#'
cutreevar <- function(obj,k=NULL,matsim=FALSE) {
  cl <- match.call()
  if (is.null(n1 <- nrow(obj$merge)) || n1 < 1) 
    stop("invalid 'tree' (merge component)")
  n <- n1 + 1
  if (is.null(k)) 
    stop("'k' must be specified")
  k <- as.integer(k)
  if (k < 1 || k > n) 
    stop(gettextf("elements of 'k' must be between 1 and %d", 
                  n), domain = NA)
  part1<-cutree(obj,k)  
  p <- length(obj$init)
  part <- rep(NA,p)
  for (i in 1:p) part[i] <- part1[obj$init[i]]  
  names(part)<-names(part1)
  res <- descript(part,obj$rec,matsim)
  res$call <- cl
  res$iter <- FALSE
  res$rec <- obj$rec
  class(res) <- "clustvar"
  return(res)
}

