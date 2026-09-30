#' @export
#' @name bootvar
#' @title  Bootstrap of individuals on a numeric matrix and on a categorical matrix
#' @description Draw a bootstrap sample from X.quanti and a bootstrap sample from X.quali
#' @param X.quanti a numeric matrix of data, or an object that can be coerced
#' to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param X.quali a categorical matrix of data, or an object that can be
#' coerced to such a matrix (such as a character vector, a factor or a data
#' frame with all factor columns).

bootvar <-
function(X.quanti=NULL,X.quali=NULL)
{
	if (!is.null(X.quanti)&& !is.null(X.quali))
	{
		n <- nrow(X.quanti)
		indice <- sample(1:n,n,replace=TRUE)
		Xboot.quanti <- X.quanti[indice,]
		Xboot.quali <- X.quali[indice,]
	}
	if (!is.null(X.quanti)&& is.null(X.quali))
	{
		n <- nrow(X.quanti)
		indice <- sample(1:n,n,replace=TRUE)
		Xboot.quanti <- X.quanti[indice,]
		Xboot.quali <- NULL
	}
	if (is.null(X.quanti)&& !is.null(X.quali))
	{
		n <- nrow(X.quali)
		indice <- sample(1:n,n,replace=TRUE)
		Xboot.quali <- X.quali[indice,]
		Xboot.quanti <- NULL
	}

	list(Xboot.quanti=Xboot.quanti,Xboot.quali=Xboot.quali)
}

