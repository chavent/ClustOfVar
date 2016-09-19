##' @export
##' @name kmeansvar
##' @title k-means clustering of variables
##' 
##' @description Iterative relocation algorithm of k-means type which performs a
##' partitionning of a set of variables.  Variables can be quantitative,
##' qualitative or a mixture of both. The center of a cluster of variables is a
##' synthetic variable but is not a 'mean' as for classical k-means. This
##' synthetic variable is the first principal component calculated by PCAmix.
##' PCAmix is defined for a mixture of qualitative and quantitative variables
##' and includes ordinary principal component analysis (PCA) and multiple
##' correspondence analysis (MCA) as special cases. The homogeneity of a
##' cluster of variables is defined as the sum of the correlation ratio (for
##' qualitative variables) and the squared correlation (for quantitative
##' variables) between the variables and the center of the cluster, which is in
##' all cases a numerical variable. Missing values are replaced by means for
##' quantitative variables and by zeros in the indicator matrix for qualitative
##' variables. 
##' 
##' @param X.quanti a numeric matrix of data, or an object that can be coerced to such a matrix 
##' (such as a numeric vector or a data frame with all numeric columns).
##' @param X.quali a categorical matrix of data, or an object that can be coerced to such a matrix 
##' (such as a character vector, a factor or a data frame with all factor columns).
##' @param init either the number of clusters or an initial partition (a vector
##' of integers indicating the cluster to which each variable is allocated).
##' If \code{init} is a number, a random set of (distinct) columns in
##' \code{X.quali} and \code{X.quanti} is chosen as the initial cluster
##' centers.
##' @param iter.max the maximum number of iterations allowed.
##' @param nstart if \code{init} is a number, \code{nstart} corresponds with
##' the number of random sets used in the process.
##' @param matsim boolean, if 'TRUE', the matrices of similarities between
##' variables in same cluster are calculated.
##' @return \item{var}{a list of matrices of squared loadings i.e. for each
##' cluster of variables, the squared loadings on first principal component of
##' PCAmix. For quantitative variables (resp. qualitative), squared loadings
##' are the squared correlations (resp. the correlation ratios) with the first
##' PC (the cluster center).} \item{sim}{a list of matrices of similarities
##' i.e. for each cluster, similarities between their variables.  The
##' similarity between two variables is defined as a square cosine: the square
##' of the Pearson correlation when the two variables are quantitative; the
##' correlation ratio when one variable is quantitative and the other one is
##' qualitative; the square of the canonical correlation between two sets of
##' dummy variables, when the two variables are qualitative. \code{sim} is
##' 'NULL' if 'matsim' is 'FALSE'.} \item{cluster}{a vector of integers
##' indicating the cluster to which each variable is allocated.} \item{wss}{the
##' within-cluster sum of squares for each cluster: the sum of the correlation
##' ratio (for qualitative variables) and the squared correlation (for
##' quantitative variables) between the variables and the center of the
##' cluster.} \item{E}{the pourcentage of homogeneity which is accounted by the
##' partition in k clusters.} \item{size}{the number of variables in each
##' cluster.} \item{scores}{a n by k numerical matrix which contains the k
##' cluster centers. The center of a cluster is a synthetic variable: the first
##' principal component calculated by PCAmix.  The k columns of \code{scores}
##' contain the scores of the n observations units on the first PCs of the k
##' clusters.} \item{coef}{a list of the coefficients of the linear
##' combinations defining the synthetic variable of each cluster.}
##' @details If the quantitative and qualitative data are in a same dataframe, the function 
##' \code{splitmix} can be used to extract automatically the qualitative and the quantitative 
##' data in two separated dataframes.
##' @author Marie Chavent \email{Marie.Chavent@@u-bordeaux.fr}, Vanessa Kuentz, Amaury Labenne,
##' Benoit Liquet, Jerome Saracco
##' @seealso
##' \code{\link{splitmix}}, \code{\link{summary.clustvar}},\code{\link{print.clustvar}},\code{\link{stability}},\code{\link{cutreevar}},\code{\link{predict.clustvar}}
##' @references Chavent, M., Liquet, B., Kuentz, V., Saracco, J. (2012),
##' ClustOfVar: An R Package for the Clustering of Variables. Journal of
##' Statistical Software, Vol. 50, pp. 1-16.
##' 
##' Chavent, M., Kuentz, V., Saracco, J. (2011), Orthogonal Rotation in PCAMIX.
##' Advances in Classification and Data Analysis, Vol. 6, pp. 131-146.
##' @keywords cluster multivariate
##' @examples
##' 
##' data(decathlon) 
##' #choice of the number of clusters
##' tree <- hclustvar(X.quanti=decathlon[,1:10])
##' stab <- stability(tree,B=60)
##' #a random set of variables is chosen as the initial cluster centers, nstart=10 times
##' part1 <- kmeansvar(X.quanti=decathlon[,1:10],init=5,nstart=10)
##' summary(part1)
##' #the partition from the hierarchical clustering is chosen as initial partition
##' part_init<-cutreevar(tree,5)$cluster
##' part2<-kmeansvar(X.quanti=decathlon[,1:10],init=part_init,matsim=TRUE)
##' summary(part2)
##' part2$sim
##' 
kmeansvar <-
  function(X.quanti=NULL, X.quali=NULL, init,iter.max=150,nstart=1,matsim=FALSE)
  {
    #init:  Either the number of clusters or a vector with group memberships
    #iter.max: The maximum number of iterations allowed.
    #nstart: If init is a number, how many random sets should be chosen?
    
    cl <- match.call()
    #split.data<-splitmix(data)
    #X.quanti<-split.data$X.quanti
    #X.quali<-split.data$X.quali
    
    rec <- recod(X.quanti,X.quali)
    
    n <- rec$n
    p <- rec$p		#total number of variables
    p1 <- rec$p1	#number of numerical variables
    X <- rec$X 		#data matrix (mixed if necessary)
    Xn <- rec$Xn	#data matrix without any replication
    Z <- rec$Z		#data matrix for input in SVD
    indexj <- rec$indexj #variables indices in columns of Z
    G <- rec$G		#dummy variables if X.quanti non null
    Gcod <- rec$Gcod	
    #	X.quanti<-rec$X.quanti
    #	X.quali<-rec$X.quali
    
    sqcc <- function(X1,X2) {
      n <- nrow(X1)
      r <- ncol(X1)
      s <- ncol(X2)
      m <- which.min(c(n,r,s))
      if (m==1) {
        A1 <- X1%*%t(X1)/n
        A2 <- X2%*%t(X2)/n
        A <- A1%*%A2
        e <- eigen(A)
        sim <- Re(e$values[2]) 
      } else {
        V12 <- t(X1)%*%X2/n
        V21 <- t(X2)%*%X1/n
        if (m==2) V<-V12%*%V21
        if (m==3) V<-V21%*%V12
        e <- eigen(V)
        sim <- Re(e$values[2]) 
      }
      return(sim)
    }	
    
    eta2 <- function(x, gpe) {
      moyennes <- tapply(x, gpe, mean)
      effectifs <- tapply(x, gpe, length)
      varinter <- (sum(effectifs * (moyennes - mean(x))^2))
      vartot <- (var(x) * (length(x) - 1))
      res <- varinter/vartot
      return(res)
    }
    
    
    partinit <- function(centers)
    {
      A <- matrix(,p,k)
      #only quantitative data
      if (p1==p) {
        A <- (t(Z)%*%Z[,centers]/n)^2
        part <- apply(A,1,which.max)
      }
      #only qualitative data
      if (p1==0) {
        for (g in 1:k) {
          X1 <- Gcod[,which(indexj==centers[g])]
          for (j in 1:p) {
            X2 <- Gcod[,which(indexj==j)]
            A[j,g] <- sqcc(X1,X2)
          }
        }
        part <- apply(A,1,which.max)
      }
      #a mixture of both
      if (p1 %in% 1:(p-1)) {
        nc1 <- length(centers[which(centers<=p1)])
        nc2 <- k-nc1
        if (nc1 != 0) {
          A[1:p1,1:nc1] <- t(Z[,1:p1])%*%Z[,centers[1:nc1]]/n
          for (g in 1:nc1) {
            for (j in (p1+1):p) A[j,g] <- eta2(X.quanti[,centers[g]],X.quali[,(j-p1)])
          }		
        }
        if (nc2 != 0) {
          for (g in (nc1+1):k) {
            for (j in 1:p1) A[j,g] <- eta2(X.quanti[,j], X.quali[,centers[g]-p1])
            X1 <- Gcod[,which(indexj==centers[g])-p1]
            for (j in (p1+1):p)  {
              X2 <- Gcod[,which(indexj==j)-p1]
              A[j,g] <- sqcc(X1,X2)
            }
          }	
        }
        part <- apply(A,1,which.max)
      }
      return(part)
    }
    
    if (missing(init)) 
      stop("'init' must be a number or a vector")
    if (length(init) == 1) {
      k <- init 
      centers <- sort(sample.int(p, k))
      part <- as.factor(partinit(centers))
    } else {
      if (!is.integer(init))
        stop("init must be a vector of integer")
      if (nstart!=1)
        stop("nstart must be equal to one")
      part <- as.factor(init)
      k <- length(levels(part))
      if (p < k) 
        stop("more cluster centers than variables")
      if (length(which(init>k))> 0)
        stop("clusters must be numbered from 1 to k")
      if (p != length(part)) 
        stop("the length of init must be equal to the number of variables")
    }
    
    do_one <- function(part)
    {
      A <- matrix(,p,k)
      iter <- 0
      diff <- TRUE
      while ((iter< iter.max) && diff)
      {
        iter <- iter+1
        
        #Representation step
        indexk<-NULL
        for (i in 1:length(indexj))
          indexk[i] <- part[indexj[i]]
        latent.var <- matrix(,n,k)
        sv<-rep(NA,k)
        for (g in 1:k) {
          Zclass <- Z[,which(indexk==g)]
          latent <- clusterscore(Zclass)
          latent.var[,g] <- latent$f
          sv[g] <- latent$sv }
        #Affectation step
        if (p1>0) {
          scorestand <- sweep(latent.var,2,STATS=sv,FUN="/")
          A[1:p1,]<-	(t(Z[,1:p1])%*%scorestand/n)^2
        }
        if (p1!=p)	{
          for (g in 1:k){
            sl.qual <- function(col) {
              eta2(latent.var[,g],col)
            }
            A[(p1+1):p,g]<-apply(X.quali,2,sl.qual)
          }	
        }
        part2 <- apply(A,1,which.max)
        
        diff <- !all(part==part2)
        part <- part2
      }
      wss <- sv^2 #within cluster sum of squares
      return(list(latent.var=latent.var,part=part,wss=wss,iter=iter))
    }
    
    res<-do_one(part)
    if (nstart >= 2) {
      best <- sum(res$wss)
      for (i in 2:nstart) {
        centers <- sort(sample.int(p, k))
        part<-as.factor(partinit(centers))
        res2 <- do_one(part)
        if ((z <- sum(res2$wss)) > best) {
          res <- res2
          best <- z}
      }
    }
    
    part <- res$part
    iter <- res$iter
    names(part) <- colnames(X)
    res <- descript(part,rec,matsim)
    res$call <- cl
    res$iter <- iter
    res$rec <- rec
    class(res) <- "clustvar"
    return(res)
  }
