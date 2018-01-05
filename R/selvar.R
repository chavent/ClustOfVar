#' @export
#' @name selvar
#' @title Selection of a given number of variables in each cluster.
#' @description This function selects in each cluster
#' a given number of variables  having the highest 
#' squared loadings. The squared loading of a variable in a cluster is its 
#' squared correlation (for numerical variable) and its correlation ratio (for categorical variable)
#' with the first PC of PCAmix applied to the variables of the cluster. 
#' @param part an  object of class \code{clustvar}
#' @param  nsel the number of variables selected in each cluster.
#' @return Returns a list where each element contains the  \code{nsel} selected variables.
#' @details If the number of variables in a cluster is smaller than \code{nsel}, 
#' all the variables of the cluster are selected
#' @examples
#' data(decathlon)
#' tree <- hclustvar(decathlon[,1:10])
#' part <- cutreevar(tree,4)
#' part$var
#' selvar(part,2) 


selvar <- function(part,nsel)
{
  if (!inherits(part, "clustvar")) 
    stop("use only with \"clustvar\" objects")
  if (nsel!=abs(as.integer(nsel)))
    strop("nsel must be a positive integer")
  sel <- function(x)
  {
    if (nrow(x) < nsel) 
      nsel <- nrow(x) 
    return(rownames(x)[1:nsel])
  }
  lapply(part$var,sel)
}

