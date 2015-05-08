## perform PCA to reveal subspecific origin of mouse samples
## <ids> = MM database IDs of samples
## <snps> = dataframe of subspecies-diagnostic markers from JPD
sso.pca <- function(ref.ids, unk.ids, snps, balance = TRUE, n.markers = 1000, ...) {
	
	require(plyr)
	require(reshape2)
	source("~/lib/util/R/db/genodb.R")
	
	if (balance)
		## sample an equal number of markers diagnostic for each subspecies; PCA plots will have a nice triangle shape
		markers <- unname(unlist( dlply(snps, .(ss), function(d) d[ sample.int(nrow(d), n.markers),"snpID" ]) ))
	else
		## use all diagnostic markers
		markers <- unname(unlist( dlply(snps, .(ss), function(d) d[ ,"snpID" ]) ))
	
	## fetch intensities for reference samples
	intens <- fetch.intensities(ref.ids, markers = markers, by = "id")
	
	## convert to matrix and do PCA
	intens.mat <- cbind( dcast(intens, sid ~ marker, value.var = "x"),
						 dcast(intens, sid ~ marker, value.var = "y") )
	X <- intens.mat[ ,-1 ]
	rownames(X) <- as.character(intens.mat[,1])
	pc <- prcomp(X, scale = TRUE)
	proj.refs <- predict(pc)
	
	## fetch intensities for query samples
	intens <- fetch.intensities(unk.ids, markers = markers, by = "id")
	
	## convert to matrix
	unk.mat <- cbind( dcast(intens, sid ~ marker, value.var = "x"),
					  dcast(intens, sid ~ marker, value.var = "y") )
	X <- unk.mat[ ,-1 ]
	rownames(X) <- as.character(unk.mat[,1])
	
	## project query samples onto PCA subspace defined by reference samples
	rez <- predict(pc, newdata = X)
	
	## return projections for both reference and query
	return( list(refs = proj.refs, samples = rez) )
	
}

## do PCA on array intensity, creating axes from <x> and (if supplied) projecting <y> onto those axes
pca.intensities <- function(x, y = NULL, formula = as.formula("sid ~ marker"), intensity.cols = c("x","y"), id.cols = c("sid"), ...) {
	
	require(reshape2)
	
	.make.intens.matrix <- function(z) {
		
		intens <- lapply(intensity.cols, function(c) {
			w <- dcast(z, formula, value.var = c)
			rownames(w) <- make.names(apply(w[,id.cols, drop = FALSE], 1, paste))
			w <- w[,colnames(w)[ !(colnames(w) %in% id.cols) ]]
			return(as.matrix(w))
		})
		intens.w <- do.call(cbind, intens)
		colnames(intens.w) <- make.names(colnames(intens.w))
		return(intens.w)
		
	}
	
	ref <- .make.intens.matrix(x)
	pc <- prcomp(ref, scale = TRUE)
	proj <- predict(pc)
	if (!is.null(y)) {
		topredict <- .make.intens.matrix(y)
		predicted <- predict(pc, newdata = topredict)
		proj <- rbind(proj, predicted)
	}
	
	attr(proj, "pca") <- pc
	attr(proj, "explained") <- pc$sdev^2/sum(pc$sdev^2)
	
	return(proj)
	
}

## perform PCA on intensities and build a tree from distance in PCA-space
pca.phylo <- function(...) {
	
	require(ape)
	proj <- pca.intensities(...)
	tree <- nj(dist(proj))
	
	return(tree)
	
}