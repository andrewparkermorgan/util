## sparse.R
## helper functions for handling large sparse matrices

## retrieve row-col indices for nonzero elements of a d[sg]CMatrix
get.sparse.indices <- function(x, ...) {
	
	if (!(inherits(x, "dsCMatrix") | inherits(x,"dgCMatrix")))
		stop("I was expecting an object of class 'dsCMatrix.'")
	
	first <- diff(x@p)
	rr <- x@i
	cc <- do.call("c", sapply(1:ncol(x), function(i) rep(i, first[i])) )
	return( cbind(rr,cc) )
	
}

## turn a sparse matrix into a dataframe of (row,col,value)
sparse.to.df <- function(x, ...) {
	
	rez <- as.data.frame(cbind(get.sparse.indices(x), x@x))
	colnames(rez) <- c("row","col","value")
	return(rez)
	
}