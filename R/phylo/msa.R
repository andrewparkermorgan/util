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

subset.alignment <- function(a, seqs, ...) {
	
	n <- 0
	aa <- list(nam = character(), nb = 0, com = NA, seq = list())
	for (s in seqs) {
		i <- match(s, a$nam)
		if (!is.na(i)) {
			aa$nam <- c(aa$nam, s)
			aa$seq <- c(aa$seq, a$seq[[i]])
			aa$nb <- aa$nb + 1
		}
	}
	class(aa) <- c("alignment", class(aa))
	return(aa)
	
}

slice <- function(a, i, ...) {
	
	newseq <- lapply(a$seq, function(x) paste(unlist(strsplit(x, ""))[ i ], collapse = ""))
	#print(newseq)
	a$seq <- newseq
	return(a)
	
}

slice2 <- function(a, ir, ...) {
	
	newseq <- lapply(a$seq, function(x) paste(unlist(strsplit(x, ""))[ (start(ir)[1]):(end(ir)[1]) ], collapse = ""))
	#print(newseq)
	a$seq <- newseq
	return(a)
	
}

.make.sdp <- function(x, ...) {
	
	code <- paste0(as.numeric(factor(x)), collapse = "")
	nalleles <- length(unique(x))
	return( list(code = code, nalleles = nalleles) )
	
}

compatible.sdps <- function(x, y, strict = FALSE, ...) {
	
	if (is.character(y))
		y <- as.numeric(strsplit(y, "")[[1]])
	
	compat <- TRUE
	if (max(y) > 1 || strict) {
		## pattern specified with 0,1,2: 0=ignore, 1=A1, 2=A2
		light <- lapply(strsplit(as.character(x), ""), function(z) unique(z[ y == 1 ]))
		dark <- lapply(strsplit(as.character(x), ""), function(z) unique(z[ y == 2 ]))
		compat <- all(sapply(light, length) == 1) && all(sapply(dark, length) == 1)
		compat <- compat && light[[1]][1] != dark[[1]][1]
	}
	else {
		light <- lapply(strsplit(as.character(x), ""), function(z) unique(z[ y == 1 ]))
		compat <- all(sapply(light, length) == 1)
	}
	
	return(compat)
	
}

beslice <- function(aln, ranges, ...) {
	
	if (!inherits(ranges, "IRanges"))
		if (is.matrix(ranges) & ncol(ranges) >= 1) {
			message("Coercing matrix input to IRanges.")
			ranges <- IRanges(start = ranges[,1], end = ranges[,2])
		}
		else
			stop("Please supply ranges as an IRanges object or a 2-column matrix.")
	
	ldply(1:length(ranges), function(i) {
		alleles <- unlist(slice2(aln, ranges[i])$seq)
		names(alleles) <- aln$nam
		rez <- data.frame(id = i, start = start(ranges)[i], end = end(ranges)[i], seq = aln$nam, allele = alleles)
		sdp <- .make.sdp(alleles)
		rez$sdp <- sdp$code
		rez$nalleles <- sdp$nalleles
		return(rez)
	})
	
}

windows <- function(a, size = 50*3, step = size, ...) {
	
	len <- length(a)
	starts <- seq(1, len, step)
	ends <- pmin(starts+size-1, len)
	return( cbind(start = starts, end = ends) )
	
}

#as.matrix.alignment <- function(a, ...) {
#	
#	seqlist <- sapply(a$seq, "[[", 1)
#	names(seqlist) <- a$nam
#	return( t(sapply(seqlist, .str.to.vector)) )
#	
#}

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

## compute omega = E[Ka/Ks] for an alignment; return matrix with rows/cols reordered by hclust
calc.omega <- function(aln, ...) {
	
	## trim alignment to multiple of 3
	trimby <- length(aln) %% 3
	if (trimby > 0)
		warning("Trimming alignment to length divisible by 3; is it really a codon alignment?")
	aln <- slice(aln, 1:(length(aln)-trimby))
	
	K <- kaks(aln)
	.kaks <- as.matrix(with(K, ka/ks))
	.kaks[ !is.finite(.kaks) | abs(.kaks) > 10 ] <- NA
	o <- hclust(dist(.kaks))$order
	.kaks <- .kaks[ o,o ]
	attr(.kaks, "order") <- o
	return(.kaks)
	
}

find.variable.ranges <- function(aln, ...) {
	
	if (inherits(aln, "alignment"))
		aln <- as.matrix.alignment(aln)
	
	var <- as.integer(unname(apply(aln, 2, function(x) length(unique(toupper(x))) > 1)))
	rl <- Rle(var)
	i <- runValue(rl) == 1
	IRanges( start = start(rl)[i], end = end(rl)[i] )
	
}