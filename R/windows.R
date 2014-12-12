make.windows <- function(gr = NULL, size = 1000, step = size, ...) {
	
	if (inherits(gr, "Seqinfo")) {
		gr <- as(gr, "GRanges")
	}
	
	windowify <- function(g) {
		
		left <- start(g)[1]
		right <- end(g)[1]
		
		starts <- seq(left, right, step)
		ends <- pmin(starts + size - 1, right)
		
		return( GRanges(seqnames = seqnames(g)[1], ranges = IRanges(start = starts, end = ends)) )
		
	}
	
	.gr <- reduce(gr, ...)
	.grl <- split(.gr, 1:length(.gr))
	rez <- endoapply(.grl, windowify)
	return( unname(unlist(rez)) )
	
}

ggply <- function(gr, ...) {
	
	require(plyr)
	
	df <- as.data.frame(gr)
	rez <- ddply(df, ...)
	rez$width <- NULL
	# print(rez)
	
	makeGRangesFromDataFrame(rez, keep.extra.columns = TRUE)
	
}

wapply <- function(gr, windows, fn, ..., unlist = TRUE) {
	
	require(plyr)
	
	olap <- as.data.frame(findOverlaps(windows, gr))
	grouped <- lapply(1:length(windows), function(w) {
		gr[ olap$subjectHits[ olap$queryHits == w ] ]
	})
	#print(length(grouped[[1]]))
	# print(class(grouped))
	# .grl <- as(grouped, "GRangesList")
	rez <- endoapply(grouped, fn, ...)
	
	if (unlist) {
		values(windows) <- cbind(values(windows), data.frame(unlist(rez)))
		return(windows)
	}
	else {
		return(rez)
	}
	
}