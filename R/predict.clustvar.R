#' @export
#' @name predict.clustvar
#' @title Scores of new objects on the synthetic variables of a given partition
#' @description A partition of variables obtained with kmeansvar or with cutreevar is given in input.
#' Each cluster of this partition is associated with a synthetic variable which is a linear combination of the variables of the cluster.
#' The coefficients of these k linear combinations (one for each cluster) are used here to calculate new scores of a objects described in a new dataset (with the same variables).
#' The output is the matrix of the scores of these new objects on the k synthetic variables.
#' @param object  an object of class clustvar
#' @param X.quanti  numeric matrix of data for the new objects
#' @param X.quali  a categorical matrix of data for the new objects
#' @param \dots Further arguments to be passed to or from other methods. They are ignored in this function.
#' @return Returns the matrix of the scores of the new objects on the k syntetic variables of the k-clusters partition given in input.
#' @examples
#' data(wine)
#' n <- nrow(wine)
#' sub <- 10:20
#' data.sub <- wine[sub,] #learning sample
#' X.quanti <- wine[sub,c(3:29)] #learning sample
#' X.quali <- wine[sub,c(1,2)]
#' part <-kmeansvar(X.quanti, X.quali, init=5)
#' X.quanti.t <- wine[-sub,c(3:29)]
#' X.quali.t <- wine[-sub,c(1,2)]
#' new <- predict(part,X.quanti.t,X.quali.t)

predict.clustvar <- function(object, X.quanti=NULL,X.quali=NULL,...)
{
  part <- object
  
  if (!inherits(part, "clustvar")) 
    stop("use only with \"clustvar\" objects")

  pfin <- part$cluster
  indexj <- part$rec$indexj
  indexg <- NULL
  
  for (i in 1:length(indexj))
    indexg[i]<-pfin[indexj[i]] 
  
  beta <- part$coef
  
 ####################### 
  train.rec <- part$rec
  train.rename.level <-  TRUE
  
  if (!is.null(X.quanti)) 
  {
    label <- rownames(X.quanti)
    n1 <- nrow(X.quanti)
    if (is.null(train.rec$X.quanti))
      stop("No quantitative dataset for training PCAmix", call. = FALSE)       
    if (!setequal(colnames(train.rec$X.quanti), colnames(X.quanti)))
      stop("The names of the columns in X.quanti and in the learning dataset
           are different", call. = FALSE)
    Y1 <- as.matrix(X.quanti[,colnames(train.rec$X.quanti)])
    # imput missing values with mean values in the train dataset
    if (any(is.na(Y1))) 
    {
      for (v in 1:ncol(Y1))
      {
        ind <- which(is.na(Y1[,v])==TRUE)
        if(length(ind)>=1)
          Y1[ind,v] <- train.rec$g[v]
      }
    }
  } else Y1 <- NULL
  
  if (!is.null(X.quali))
  {
    label <- rownames(X.quali)
    n2 <- nrow(X.quali)
    p2 <- ncol(X.quali)
    if (is.null(train.rec$X.quali))
      stop("No qualitative dataset for training PCAmix", call. = FALSE)  
    if (!setequal(colnames(train.rec$X.quali), colnames(X.quali)))
      stop("The names of the columns in X.quali and in the learning dataset
           are different", call. = FALSE)
    GNA <- PCAmixdata::tab.disjonctif.NA(X.quali, 
                                         rename.level = train.rename.level)
    G <- as.matrix(replace(GNA,is.na(GNA),0))

    if (!setequal(colnames(train.rec$G),colnames(G)))
    {
      #levels in train not in test : levels are merged in test
      if (length(setdiff(colnames(train.rec$G), colnames(G))) > 0)
      {
        for (v in 1:p2)
          levels(X.quali[,v]) <- union(levels(X.quali[,v]),
                                              levels(train.rec$X.quali[,v]))
        GNA <- PCAmixdata::tab.disjonctif.NA(X.quali, 
                                             rename.level = train.rename.level)
        G <- as.matrix(replace(GNA,is.na(GNA),0))
      }
      
      #levels in test not in train : observations are removed
      if (length(setdiff(colnames(G), colnames(train.rec$G))) > 0)
      { 
        if (length(setdiff(colnames(G), colnames(train.rec$G))) == length(colnames(G)))
          stop("No level in common between the test and the training dataset",
               call. = FALSE)
        G2 <- G[,is.element(colnames(G), colnames(train.rec$G)), drop=FALSE]
        if (length(which(apply(G2,1,sum) < p2)) == nrow(G))
          stop("No observation in the test dataset can be predicted", 
               call. = FALSE)
        G <- G2
        if (any(apply(G2,1,sum) < p2))
        {
          warning("Predictions can not be preformed for some observations")
          G <- G2[-which(apply(G2,1,sum) < p2),, drop=FALSE]
          if (!is.null(Y1))
            Y1 <- Y1[-which(apply(G2,1,sum) < p2),,drop=FALSE]
        }
      }
    }
  G <- G[ , colnames(train.rec$G), drop = FALSE]
  } else G <- NULL
  
  if (!is.null(X.quanti) && !is.null(X.quali))
  {
    if (n1 != n2) 
      stop("The number of objects in X.quanti and X.quali must be 
                       the same", call. = FALSE)
    if (sum(rownames(X.quali)!=rownames(X.quanti))!=0) 
      stop("The names/order of the objects in X.quanti and X.quali must be the same", 
           call. = FALSE)
  }
  
  Y <- cbind(Y1,G)
  n <- nrow(Y)
  
  scores <- matrix(NA, nrow(Y), length(beta))
  
  for (g in 1: length(beta))
  {
    Yg <- Y[,which(indexg==g), drop=FALSE]
    scores[,g] <-Yg %*% beta[[g]][-1] +  beta[[g]][1]
  }
  colnames(scores) <- paste("cluster", 1:length(beta), sep = "")
  rownames(scores) <- rownames(Y)
  return(scores)
}
