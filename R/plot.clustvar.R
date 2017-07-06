#' @export plot.clustvar
#' @export
#' @name plot.clustvar
#' @title  Plot loadings in each cluster.
#' @description Plot dotchart with the "loadings" of the variables in each cluster. The loading of a
#' numerical variable is the correlation between this variables and the synthetic variable of its cluster.
#' The loading of the level of a categorical variable is the mean value of the synthetic variable
#' of the cluster on observations having this level.
#' @param x an object of class \code{clustvar} obtained with \code{cutreevar} or \code{kmeansvar}.
#' @param \dots Further arguments to be passed to or from other methods. They
#' are ignored in this function.
#' @return \item{coord.quanti}{coordinates of quantitative variables belonging to cluster k on the synthetic
#'  variable associate to the same cluster k}
#' \item{coord.levels}{coordinates of levels of categorical variables belonging to cluster k on the synthetic
#'  variable associate to the same cluster k}
#' @examples
#' data(wine)
#' X.quanti <- PCAmixdata::splitmix(wine)$X.quanti
#' X.quali <- PCAmixdata::splitmix(wine)$X.quali
#' tree <- hclustvar(X.quanti,X.quali)
#' tree.cut<-cutreevar(tree,6)
#'
#' #plot of scores on synthetic variables
#' res.plot <- plot(tree.cut)
#' res.plot$coord.quanti
#' res.plot$coord.levels

plot.clustvar<-function(x,...){
  
  obj<-x
  base.tot <- obj$rec$X
  parti<-obj$cluster
  names.group <- names(obj$var)
  k<-length(names.group)
  
  base.cluster<-PCAmixdata::splitgroups(data=base.tot, groups=parti, name.groups=names.group)$data.groups
  res.coord.quanti<-list()
  for (i in 1:k){res.coord.quanti[[i]]<-"No quantitative variables in this cluster"}
  res.coord.levels<-list()
  for (i in 1:k){res.coord.levels[[i]]<-"No categorical variables in this cluster"}
  
  for (i in 1:k){
    cutmix<-PCAmixdata::splitmix(base.cluster[[i]])
    res.PCAmix<-PCAmixdata::PCAmix(X.quanti=cutmix$X.quanti, X.quali=cutmix$X.quali, rename.level=TRUE, graph=FALSE)
    typ.group="MIX"
    if (is.null(cutmix$X.quali)) 
      typ.group="QT"
    if (is.null(cutmix$X.quanti)) 
      typ.group="QL"
    if (typ.group=="QT"){
      coord.quanti <- res.PCAmix$quanti$coord[, 1,drop=F]
      coord.quanti <- coord.quanti[order(coord.quanti), ,drop=F]
      dotchart(coord.quanti[, 1, drop=F],xlab="Correlation with the synthetic variable", 
               main=paste0("Cluster ",i), xlim=c(-1,1),gcolor="white")
      res.coord.quanti[[i]] <- coord.quanti[, 1, drop=F]
      #plot(res.PCAmix, choice="cor", axes=c(1,2),main=paste0("Correlation circle of quantitative variables of ",names.group[i]))
    }
    
    if (typ.group=="QL"){
      coord.categ <- res.PCAmix$levels$coord[, 1,drop=F]
      coord.categ <- coord.categ[order(coord.categ), ,drop=F]      
      dotchart(coord.categ[, 1, drop=F],xlab="Mean value of the synthetic variable", 
               main=paste0("Cluster ",i),gcolor="white")
      res.coord.levels[[i]] <- coord.categ[, 1, drop=F]
    }
    
    if (typ.group=="MIX"){
      coord.categ <- res.PCAmix$levels$coord[, 1,drop=F]
      coord.categ <- coord.categ[order(coord.categ), ,drop=F]      
      dotchart(coord.categ[, 1, drop=F],xlab="Mean value of the synthetic variable", 
               main=paste0("Cluster ",i),gcolor="white")
      res.coord.levels[[i]] <- coord.categ[, 1, drop=F]
      
      
      coord.quanti <- res.PCAmix$quanti$coord[, 1,drop=F]
      coord.quanti <- coord.quanti[order(coord.quanti), ,drop=F]      
      dotchart(coord.quanti[, 1, drop=F],xlab="Correlation with the synthetic variable", 
               main=paste0("Cluster ",i), xlim=c(-1,1),gcolor="white")  
      res.coord.quanti[[i]] <- coord.quanti[, 1, drop=F]
      #plot(res.PCAmix, choice="cor", axes=c(1,2),cex=0.5, main=paste0("Correlation circle of quantitative variables of ",names.group[i]))
    }
  }
  names(res.coord.quanti) <- names (res.coord.levels) <- paste("Cluster",1:k,sep="")
  res <- list(coord.quanti=res.coord.quanti, coord.levels=res.coord.levels)
  return(res)
}
