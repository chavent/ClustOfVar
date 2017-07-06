#' @export
#' @name descript
#' @title Description of a partition of variables
#' @description Gives a description of a partition of variables and the scores of the cluster's synthetic variables.
#' @param part a vector with cluster memberships of a partition of variables
#' @param rec the element $rec of an object of class clustvar
#' @param  matsim boolean, if TRUE, the matrices of similarities between variables in same cluster are calculated.
#' @return 
#' \item{var}{a list of matrices of squared loadings i.e. for each cluster of variables, the squared loadings on first principal component of PCAmix.
#' For quantitative variables (resp. qualitative), squared loadings are the squared correlations (resp. the correlation ratios) with the first PC (the cluster center).}
#' \item{sim}{a list of matrices of similarities i.e. for each cluster, similarities between their variables. 
#' The similarity between two variables is defined as a square cosine: the square of the Pearson correlation when the two variables are quantitative; 
#' the correlation ratio when one variable is quantitative and the other one is qualitative; 
#' the square of the canonical correlation between two sets of dummy variables, when the two variables are qualitative.
#' \code{sim} is 'NULL' if 'matsim' is 'FALSE'.}
#' \item{cluster}{a vector of integers indicating the cluster to which each variable is allocated.}
#' \item{wss}{the within-cluster sum of squares for each cluster: the sum of the correlation ratio (for qualitative variables) and the squared correlation (for quantitative variables) between the variables and the center of the cluster.}
#' \item{E}{the pourcentage of  homogeneity which is accounted by the partition in k clusters.}
#' \item{size}{the number of variables in each cluster.}
#' \item{scores}{a n by k numerical matrix which contains the k cluster centers. The center of a cluster is a synthetic variable: the first principal component calculated by PCAmix. 
#' The k columns of \code{scores} contain the scores of the n observations units on the first PCs of the k clusters.}
#' \item{coef}{a list of the coefficients of the linear combinations defining the synthetic variable of each cluster.}
#' @keywords internal
descript<-function(part,rec,matsim=FALSE) 
{  
  k<-length(levels(as.factor(part)))
  
  Z<-rec$Z
  X<-rec$X  
  gc <- rec$g 
  sdev <- rec$s 
  p1<-rec$p1	
  indexj<-rec$indexj 
  n<-rec$n
  p<-rec$p
  if (p!=length(part)) 
    stop("The data matrix recoded in 'rec' must have the same number of variables than the partition")
  indexk<-NULL
  for (i in 1:length(indexj))
    indexk[i]<-part[indexj[i]]
  latent.var <- matrix(,n,k)
  sv<-rep(NA,k)#singular values
  
  coef <- structure(vector(mode = "list", length = k), names = paste("cluster", 1:k, sep = ""))
  var <- structure(vector(mode = "list", length = k), names = paste("cluster", 1:k, sep = ""))
  sim <- structure(vector(mode = "list", length = k), names = paste("cluster", 1:k, sep = ""))
  
  for (g in 1:k)  {
    Zclass<-as.matrix(Z[,which(indexk==g)])
    gclass <- gc[which(indexk==g)]
    sclass <- sdev[which(indexk==g)]
    latent<-clusterscore(Zclass)
    latent.var[,g]<- latent$f
    sv[g]<-latent$sv
    v <- 	latent$v
    clus<-which(part==g)
    #correlation ratio
    C <- matrix(NA, nrow=length(clus), ncol=2)
    colnames(C)<-c("squared loading","correlation")
    rownames(C)<-names(clus)
    for (i in 1:length(clus)){
      C[i,1] <- mixedVarSim(latent.var[,g],X[,clus[i]])
      if (is.numeric(X[,clus[i]])) #ajout
        C[i,2]<-cor(X[,clus[i]], latent.var[,g], use="complete.obs")#ajout
    }
    #C<-C[order(abs(C[,2]),decreasing=TRUE), ] #ajout
    C<-C[order(C[,1],decreasing=TRUE), ] #ajout
    var[[g]] <- C
    #coefficients of the score function
    beta <- matrix(NA,length(v)+1,1)
    beta[1,1] <- -sum(v*gclass/sclass)
    beta[2:(length(v)+1),1] <- v/sclass
    rownames(beta) <- c("const",colnames(Z)[which(indexk==g)])
    coef[[g]] <- beta
    #similarity matrix
    tabcos2 <- matrix(NA, length(clus), length(clus))
    colnames(tabcos2)<-names(clus)
    rownames(tabcos2)<-names(clus)
    diag(tabcos2)<-1
    if ((length(clus)>1) && (matsim==TRUE)) {
      for (i in 1:(length(clus)-1))
        for(j in (i+1):length(clus))  {
          tabcos2[i,j] <- mixedVarSim(X[,clus[i]],X[,clus[j]])
          tabcos2[j,i]<- tabcos2[i,j]
        }
    }		
    if (matsim==TRUE) sim[[g]] <- tabcos2 else sim[[g]] <- NULL
  }
  
  wss<-sv^2 
  size<-rep(NA,k)
  for (g in 1:k) size[g]<-length(which(part==g))
  colnames(latent.var) <- names(wss) <- names(size) <- paste("cluster", 1:k, sep = "")
  rownames(latent.var) <-rownames(X)
  Ht<-clusterscore(Z)$sv^2
  E<-(sum(wss)-Ht)/(p-Ht)*100
  return(list(var=var,coef=coef,sim=sim,cluster=part,wss=wss,E=E,size=size,scores=latent.var,rec=rec,k=k))
}
