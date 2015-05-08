library(GenomicRanges)

# ## project one set of intervals onto another
# project.onto <- function(x, onto, keep = c("biggest","first","last","arbitrary"), resolve = TRUE, ...) {
# 	
# 	if (!isDisjoint(onto))
# 		stop("Target interval(s) not disjoint; projection is ambiguous.")
# 	
# 	if (!isDisjoint(x))
# 		if(!resolve)
# 			stop("Intervals to project are not disjoint.")
# 		else
# 			x <- resolve.overlaps(x)
# 	
# 	proj <- onto
# 	keep <- match.arg(keep)
# 	olap <- findOverlaps(x, onto, type = "any")
# 	sh <- subjectHits(olap)
# 	qh <- queryHits(olap)
# 	
# 	## check if any source interval projects to >1 target interval
# 	tbl <- table(subjectHits(olap))
# 	if (any(tbl > 1)) {
# 		warning(paste0("At least one target interval spans multiple source intervals; resolving using rule '", keep,"'..."))
# 		if (keep == "biggest") {
# 			## keep the metadata from the longest source interval
# 			j <- as.numeric(names(tbl[ tbl > 1 ]))
# 			tofix <- subjectHits(olap) %in% j
# 			ssh <- sh[tofix]
# 			old.sh <- sh
# 			old.qh <- qh
# 			sh <- numeric(0)
# 			qh <- numeric(0)
# 			for (i in sort(unique(ssh))) {
# 				sh <- c(sh, i)
# 				these.x <- old.qh[ old.sh == i ]
# 				biggest <- which.max(width(x[these.x]))
# 				qh <- c(qh, these.x[biggest])	
# 			}
# 			sh <- c(sh, old.sh[!tofix])
# 			qh <- c(qh, old.qh[!tofix])
# 		}
# 		else {
# 			## use GRanges built-ins to get max. 1 overlap per source interval
# 			sh <- findOverlaps(x, onto, select = keep)
# 			qh <- seq_along(x)
# 		}
# 	}
# 		
# 	## add metadata using rules defined above
# 	for (col in colnames(values(x))) {
# 		values(proj)[,col] <- NA
# 		values(proj)[ sh,col ] <- as.character(values(x)[ qh,col ])
# 		if (is.factor(values(x[,col])))
# 			values(proj)[,col] <- factor(values(proj)[,col], levels = levels(values(x)[,col]))
# 	}
# 	return(proj)
# 	
# }

## project one set of intervals onto the other, via the union of endpionts of both sets
project.onto <- function(x, onto, resolve = FALSE, seal = FALSE, ...) {
	
	if (!isDisjoint(onto))
		stop("Target interval(s) not disjoint; projection is ambiguous.")
	
	if (seal)
		## assume query and target are supposed to cover all of genome
		x <- seal.gaps(x)
	
	if (resolve)
		## force ranges in query to be non-overlapping
		x <- resolve.overlaps(x)
	
	pieces <- disjoin( c(onto, x, ignore.mcols = TRUE), ignore.strand = TRUE )
	pieces <- pieces[ pieces %over% onto & pieces %over% x ]
	cols <- colnames(values(onto))
	values(pieces) <- mergeByOverlaps(pieces, onto)[ ,cols ]
	
	olap <- mergeByOverlaps(pieces, x)
	proj <- pieces
	values(proj) <- olap[ ,cols ]
	
	return(proj)
	
}

## split a GRanges into a GRangesList by the value of specified column
split.by.column <- function(gr, col, ...) {
	
	if (!(col %in% colnames(values(gr))))
		stop("Can't find that column.")
	
	rez <- split(gr, values(gr)[ ,col ], drop = FALSE)
	
	return(rez)
	
}

## map rows in <gr> to elements of <onto> by value of <col>, then project
##	rows of <gr> onto the rows in that element
map.and.project <- function(gr, onto, col = "strain", ...) {
	
	grl <- split.by.column(gr, col)
	if (!inherits(onto, "GRangesList"))
		onto <- split.by.column(onto, col = col)
	
	mapply(project.onto, grl, onto, ...)
	
}

## resolve overlapping intervals into disjoint intervals
resolve.overlaps <- function(gr) {
	
	#print(gr)
	gr <- sort(gr)
	olap <- findOverlaps(gr, ignoreSelf = TRUE, ignoreRedundant = TRUE)
	ii <- unique(queryHits(olap))
	for (i in ii) {
		matches <- as.data.frame(olap[ queryHits(olap) == i ])
		if (nrow(matches) == 1) {
			ends <- c(end(gr[ matches[1,1] ]), start(gr[ matches[1,2] ]))
			midp <- floor(mean(range(ends)))
			end(gr)[ matches[1,1] ] <- midp
			start(gr)[ matches[1,2] ] <- midp+1
		}
		else {
			print(matches)
			stop("Complicated overlap; I probably cant solve it.")
		}
	}
	
	return(gr)
	
}

## make intervals cover whole chromosomes, doing our best to resolve overlaps
seal.gaps <- function(.gr, ...) {
	
	rez <- GRanges(seqinfo = seqinfo(.gr))
	for (c in seqlevels(.gr)) {
		
		sealed <- tryCatch({
			gr <- sort(.gr[ as.vector(seqnames(.gr)) == c ])
			if (!isDisjoint(gr)) {
				warning("Warning: ranges are not disjoint; 'sealing' may have no meaning.")
				gr <- resolve.overlaps(gr)
				#print(gr)
			}
			start(gr)[1] <- 1
			end(gr)[ length(gr) ] <- seqlengths(gr)[c]
			if (length(gr) > 1) {
				starts <- start(gr)[-1]
				ends <- end(gr)[ -length(gr) ]
				newends <- floor((starts+ends)/2)
				newstarts <- newends+1
				end(gr)[ -length(gr) ] <- newends
				start(gr)[ -1 ] <- newstarts
			}
			gr
		},
		error = function(e) {
			cat("\tsomething went wrong; moving on\n")
			GRanges(seqinfo = seqinfo(.gr))
		})
		
		rez <- c(rez, sealed)
	}
	return(trim(rez))
	
}

make.windows <- function(gr = NULL, size = 1000, step = size, ...) {
	
	si <- gr
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
	rez <- unname(unlist(rez))
	if (inherits(si, "Seqinfo"))
		seqinfo(rez) <- si
	return(rez)
	
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