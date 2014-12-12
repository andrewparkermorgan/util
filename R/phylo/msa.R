## helper functions for handling multiple sequence alignments in R

library(seqinr) # provides "alignment" class
library(GenomicRanges) # useful helper functions for run-length encoding

.str.to.vector <- function(x, split = "", ...) {
	unlist(strsplit(x, split, ...))
}

length.alignment <- function(a, ...) {
	lens <- pmax( sapply(a$seq, nchar), sapply(a$seq, length) )
	if (!all(lens == max(lens)))
		warning("Sequences not all of equal length; alignment object is messed up.")
	return(max(lens))
}

slice <- function(a, i, ...) {
	
	newseq <- lapply(a$seq, function(x) paste(unlist(strsplit(x, ""))[ i ], collapse = ""))
	#print(newseq)
	a$seq <- newseq
	return(a)
	
}

windows <- function(a, size = 50*3, step = size, ...) {
	
	len <- length(a)
	starts <- seq(1, len, step)
	ends <- pmin(starts+size-1, len)
	return( cbind(start = starts, end = ends) )
	
}

as.matrix.alignment <- function(a, ...) {
	
	seqlist <- sapply(a$seq, "[[", 1)
	names(seqlist) <- a$nam
	return( t(sapply(seqlist, .str.to.vector)) )
	
}

pct.gaps <- function(a, ...) {
	amat <- as.matrix(a)
	apply(amat, 2, function(x) sum(x == "-")/length(x))
}

gaps <- function(a, threshold = 0, ...) {
	which(pct.gaps(a) > threshold)
}

ungap <- function(x, gap.char = c("-"), ...) {
	x[ !(x %in% gap.char) ]
}

## project from coordinates of a constituent sequence, onto the coordinates of the alignment
## <seq> = a possibly-gapped sequence from the alignment
## <coords> = coordinates to project
project <- function(seq, coords, gap.char = c("-"), ...) {
	
	if (length(seq) == 1 & nchar(seq[1]) > 0)
		seq <- unlist(strsplit(seq, ""))
	
	## use run-length encoding to get coordinates of gapless sequence
	seq[ !(seq %in% gap.char) ] <- "X"
	runs <- Rle(seq)
	i <- which((runValue(runs) == "X"))
	starts <- start(runs)[i]
	ends <- end(runs)[i]
	widths <- width(runs)[i]
	
	starts.adj <- rep(0, length(starts))
	ends.adj <- rep(0, length(ends))
	z <- 0
	for (j in 1:(length(starts))) {
		z.new <- z+widths[j]
		starts.adj[j] <- z + 1
		ends.adj[j] <- z.new
		z <- z.new
	}
	
	#print(cbind(starts, ends, starts.adj, ends.adj, widths, cumsum(widths)))
	ivls <- c(starts.adj, max(ends.adj))
	where <- findInterval(coords, ivls, rightmost.closed = FALSE, all.inside = TRUE)
	offset <- coords - cumsum(widths)[where-1]
	coords.adj <- starts[where] + offset - 1
	
	return(coords.adj)
	
}