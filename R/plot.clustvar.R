##' @export
##' @name plot.clustvar
##' @method plot clustvar
##' @title  Plot scores of variables on synthetic variables.
##' @description Plot scores of variables of cluster k on the synthetic variable associated to the same cluster. For each cluster
##' two plots are displayed: one for numerical variables (if available in the cluster) and one for levels of categorical
##' variables (if available in the cluster).
##' @param x an object of class \code{clustvar} obtained with \code{cutreevar} or \code{kmeansvar}.
##' @param \dots Further arguments to be passed to or from other methods. They
##' are ignored in this function.
##' @return \item{coord.quanti}{coordinates of quantitative variables belonging to cluster k on the synthetic
##'  variable associate to the same cluster k}
##' \item{coord.levels}{coordinates of levels of categorical variables belonging to cluster k on the synthetic
##'  variable associate to the same cluster k}
##' @keywords cluster multivariate
##' @examples
##' 
##' data(wine)
##' X.quanti <- splitmix(wine)$X.quanti
##' X.quali <- splitmix(wine)$X.quali
##' tree <- hclustvar(X.quanti,X.quali)
##' tree.cut<-cutreevar(tree,6)
##' 
##' #plot of scores on synthetic variables
##' res.plot<-plot(tree.cut) 
##' res.plot$coord.quanti
##' res.plot$coord.levels
plot.clustvar<-function(x,...){
  
  obj<-x
  base.tot <- obj$rec$X
  parti<-obj$cluster
  names.group<-names(obj$var)
  k<-length(names.group)
  
  base.cluster<-splitgroups(base=base.tot, groups=parti, name.groups=names.group)
  res.coord.quanti<-list()
  for (i in 1:k){res.coord.quanti[[i]]<-"No quantitative variables in this cluster"}
  res.coord.levels<-list()
  for (i in 1:k){res.coord.levels[[i]]<-"No categorical variables in this cluster"}
  
  for (i in 1:k){
    cutmix<-splitmix(base.cluster[[i]])
    res.PCAmix<-PCAmix(X.quanti=cutmix$X.quanti, X.quali=cutmix$X.quali, rename.level=TRUE, graph=FALSE)
    typ.group<-cutmix$'typ.group'
    
    if (typ.group=="QT"){
      coord.quanti <- res.PCAmix$quanti$coord[, 1,drop=F]
      coord.quanti <- coord.quanti[order(coord.quanti), ,drop=F]
      dotchart(coord.quanti[, 1, drop=F],xlab=paste0("Synthetic variable of cluster ",i), 
               main="Correlations between quantitative variables and the synthetic variable", xlim=c(-1,1),gcolor="white",...)
      res.coord.quanti[[i]] <- coord.quanti[, 1, drop=F]
      #plot(res.PCAmix, choice="cor", axes=c(1,2),main=paste0("Correlation circle of quantitative variables of ",names.group[i]))
    }
    
    if (typ.group=="QL"){
      coord.categ <- res.PCAmix$levels$coord[, 1,drop=F]
      coord.categ <- coord.categ[order(coord.categ), ,drop=F]      
      dotchart(coord.categ[, 1, drop=F],xlab=paste0("Synthetic variable of cluster ",i), 
               main="Coordinates of levels of qualitative variables on the synthetic variable",gcolor="white",...)
      res.coord.levels[[i]] <- coord.categ[, 1, drop=F]
    }
    
    if (typ.group=="MIX"){
      coord.categ <-res.PCAmix$levels$coord[, 1,drop=F]
      coord.categ <- coord.categ[order(coord.categ), ,drop=F]      
      dotchart(coord.categ[, 1, drop=F],xlab=paste0("Synthetic variable of cluster ",i), 
               main="Coordinates of levels of qualitative variable on the synthetic variable",gcolor="white",...)
      res.coord.levels[[i]] <- coord.categ[, 1, drop=F]
      
      
      coord.quanti <- res.PCAmix$quanti$coord[, 1,drop=F]
      coord.quanti <- coord.quanti[order(coord.quanti), ,drop=F]      
      dotchart(coord.quanti[, 1, drop=F],xlab=paste0("Synthetic variable of cluster ",i), 
               main="Correlations between quantitative variables and the synthetic variable", xlim=c(-1,1),gcolor="white",...)  
      res.coord.quanti[[i]] <- coord.quanti[, 1, drop=F]
      #plot(res.PCAmix, choice="cor", axes=c(1,2),cex=0.5, main=paste0("Correlation circle of quantitative variables of ",names.group[i]))
    }
  }
  names(res.coord.quanti) <- names (res.coord.levels) <- paste("Cluster",1:k,sep="")
  res<-list(coord.quanti=res.coord.quanti, coord.levels=res.coord.levels)
  return(res)
}
