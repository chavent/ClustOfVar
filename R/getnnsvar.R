##' @export
##' @name getnnsvar
##' @title Nearest neighbor of variables 
##' @description Nearest neighbor of variables 
##' @param diss a dissimilarity matrix between variables
##' @param flag a vector of size \code{p} which indicates if we want to compute nearest
##'  neighbor of variable j (flag[j]=1) or not (flag[j]=0)
##' @keywords internal


getnnsvar <-
  function(diss, flag) {
    # Inputs.  diss: full distance matrix.
    #          flag: "live" rows indicated by 1 are to be processed.
    # Returns. List of: nn, nndiss.
    #          nn:   list of nearest neighbor of each row.
    #          nndiss: nearest neigbbor distance of each row.
    
    mini.vec.flag<-function(vect_flag){
      L<-length(vect_flag)
      flagi<-vect_flag[L]
      vect<-vect_flag[-L]
      if (flagi==0){res<-c(0,0)}
      else{
        search.min<-which(vect==min(vect[vect>0]))[1]
        minobs<-search.min
        mindis<-vect[search.min]
        res<-c(round(minobs,1),mindis)
      }
      return(res)                          
    }
    
    res<-apply(cbind(diss,flag),1,FUN=mini.vec.flag)
    result<-list(nn=res[1,],nndiss=res[2,])
    return(result)
  }



