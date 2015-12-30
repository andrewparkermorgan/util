## utility functions for handling genotype matrices

## convert between styles of chromosome names (I prefer UCSC)
convert.names <- function(chrs, to = c("ncbi","ensembl","ucsc","plink"), muga.style = FALSE, ...) {

	ucsc <- paste0("chr", c(1:19, "X","Y","M"))
	ncbi <- c(1:19,"X","Y","MT")

	to <- match.arg(to)
	converter <- function(x) { identity(x) }
	if (to == "ncbi" | to == "ensembl") {
		converter <- function(x) setNames(ncbi, ucsc)[x]
	}
	else if (to == "ucsc") {
		converter <- function(x) setNames(ucsc, ncbi)[x]
	}
	else if (to == "plink") {
		#message("converting to plink")
		converter <- function(x) {
			nn <- gsub("^chr","", x)
			nn <- gsub("^M$","MT", nn)
			return(nn)
		}
	}

	nn <- converter( as.character(chrs) )
	if (muga.style)
		nn <- gsub("^MT$","M", chrs)

	return(nn)

}

clean.names <- function(x, space = "", ...) {

	x <- gsub("^\\s+|\\s+$", "", x)
	x <- make.unique(x)
	x <- gsub(" ", space, x)
	return(x)

}

geno.to.matrix <- function(gty, map = NULL, use.names = "id", by = "marker", make.names = FALSE, value.col = "call", ...) {

	require(reshape2)

	## if map specified, attach it to get cM positions
	if (.is.valid.map(map)) {
		rownames(map) <- as.character(map$marker)
		if ("cM" %in% colnames(gty))
			gty[ ,cM := NULL ]
		if ("chr" %in% colnames(gty))
			gty[ ,chr := NULL ]
		if ("pos" %in% colnames(gty))
			gty[ ,pos := NULL ]
		if (make.names)
			gty$marker <- make.names(gty$marker)
		before <- unique(gty$marker)
		nsnps.before <- length(unique(gty$marker))
		.map <- data.table(map[ ,c("chr","marker","cM","pos") ])
		setkey(.map, "marker")
		setkey(gty, "marker")
		message(paste("Attaching map: started with", nsnps.before, "markers"))
		gty <- gty[.map]
		nsnps.after <- length(unique(gty$marker))
		message(paste("Done attaching map: ended with", nsnps.after, "markers"))
		if (nsnps.before != nsnps.after) {
			message(paste("Dropped the following markers:"))
			message( paste(setdiff(before, gty$marker) , collapse = ",") )
		}
	}
	#print(head(gty))

	## if genetic position included, keep it
	fm <- "chr + marker + pos ~ sid"
	cols <- c("chr","marker","pos")
	if ("cM" %in% colnames(gty)) {
		fm <- "chr + marker + cM + pos ~ sid"
		cols <- c("chr","marker","cM","pos")
	}
	#map.cols <- c("chr","marker","cM","pos")

	setkeyv(gty, c("sid", use.names))
	ind <- unique(gty)
	#print(head(ind))
	ind.names <- setNames( as.character(ind[ ,eval(parse(text = use.names)) ]), as.character(ind$sid) )

	gty.mat <- dcast(gty, as.formula(fm), value.var = value.col)
	gty.mat <- gty.mat[ with(gty.mat, order(chr, pos)), ]
	rownames(gty.mat) <- gty.mat$marker
	colnames(gty.mat) <- c(cols, ind.names[ colnames(gty.mat)[ (length(cols)+1):ncol(gty.mat) ] ])
	colnames(gty.mat) <- gsub(" ","", colnames(gty.mat))
	
	map <- gty.mat[ ,cols ]
	
	gty.mat <- as.matrix(gty.mat[ ,-(1:length(cols)) ])
	class(gty.mat) <- c("genotypes", class(gty.mat))
	attr(gty.mat, "map") <- map
	return(gty.mat)

}

## do some operator overloading

`[.genotypes` <- function(x, i = TRUE, j = TRUE, drop = FALSE, ...) {
	
	r <- NextMethod("[")
	if (!is.null(attr(x, "map"))) {
		attr(r, "map") <- attr(x, "map")[ i, ]
	}
	if (!is.null(attr(x, "ped"))) {
		attr(r, "ped") <- attr(x, "ped")[ j, ]
		rownames(attr(r, "ped")) <- as.character(attr(r, "ped")$iid)
	}
	if (.has.valid.intensity(x)) {
		x.new <- attr(x, "intensity")$x[ i,j, drop = FALSE ]
		y.new <- attr(x, "intensity")$y[ i,j, drop = FALSE ]
		attr(r, "intensity") <- list(x = x.new, y = y.new)
		if (!is.null(attr(x, "normalized")))
			attr(r, "normalized") <- attr(x, "normalized")
	}

	if (!is.null(attr(x, "alleles")))
		attr(r, "alleles") <- attr(x, "alleles")
	
	if (!is.null(attr(x, "filter.sites")))
		attr(r, "filter.sites") <- attr(x, "filter.sites")[i]
	if (!is.null(attr(x, "filter.samples")))
		attr(r, "filter.samples") <- attr(x, "filter.samples")[j]
	
	class(r) <- c("genotypes", class(r))
	return(r)
}

`$.genotypes` <- function(x, expr, ...) {
	
	attributes(x)[[expr]]
	
}


## set some S3 generics for useful accessor functions

markers <- function(x) UseMethod("markers")
markers.genotypes <- function(gty, ...) {
	attr(gty, "map")
}

samples <- function(x) UseMethod("samples")
samples.genotypes <- function(gty, ...) {
	if (!is.null(attr(gty, "ped")))
		attr(gty, "ped")
	else
		colnames(gty)
}

filters <- function(x) UseMethod("filters")
filters.genotypes <- function(gty, ...) {
	get.filters(gty, ...)
}

intensity <- function(x) UseMethod("intensity")
intensity.genotypes <- function(gty, ...) {
	attr(gty, "intensity")
}

## add a method to overload subset() on a genotypes object which
## preserves attributes (markers, samples, intensity...)
subset.genotypes <- function(gty, expr, by = c("markers","samples"), ...) {

	if (!inherits(gty, "genotypes"))
		stop("Only willing to subset an object of class 'genotypes.'")

	by <- match.arg(by)
	e <- substitute(expr)
	if (by == "markers") {
		if (!.has.valid.map(gty))
			stop("Can't subset genotypes by marker without a marker map.")
		r <- eval(e, attr(gty, "map"), parent.frame())
		r <- r & !is.na(r)
		return( gty[ r, ] )
	}
	else if (by == "samples") {
		if (is.null(attr(gty, "ped")))
			stop("Can't subset genotypes by samples without sample information.")
		r <- eval(e, attr(gty, "ped"), parent.frame())
		r <- r & !is.na(r)
		return( gty[ ,r ] )
	}
	else {
		stop()
	}

}

## overload cbind() to join not only genotype matrix but also sample metadata
## NB: intensities and filters are dropped at this point
cbind.genotypes <- function(a, b, ...) {
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to bind two objects of class 'genotypes.'")
	
	if (nrow(a) != nrow(b) | any(rownames(a) != rownames(b)))
		stop("Number and names of markers don't match.  Try merging genotypes instead of cbind-ing.")
	
	message(paste0("Adding ",ncol(b)," individuals to the existing ",ncol(a),"."))
	
	rez <- cbind(unclass(a), unclass(b))
	class(rez) <- c("genotypes", class(rez))
	if (!is.null(attr(a, "map")))
		attr(rez, "map") <- attr(a, "map")
	if (!is.null(attr(a, "ped")))
		attr(rez, "ped") <- rbind( attr(a, "ped"), attr(b, "ped") )
	
	return(rez)
	
}

## overload rbind() to join not only genotype matrix but also marker metadata
## NB: intensities and filters are dropped at this point
rbind.genotypes <- function(a, b, ...) {
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to bind two objects of class 'genotypes.'")
	
	cols.a <- colnames(a)
	cols.b <- colnames(b)
	if (length(setdiff(a,b)) || length(setdiff(b,a)))
		stop("Number and names of markers don't match.  Try merging genotypes instead of cbind-ing.")
	
	message(paste0("Adding ",nrow(b)," markers to the existing ",nrow(a),"."))
	
	rez <- rbind(unclass(a)[ ,cols.a ], unclass(b)[ ,cols.a ])
	class(rez) <- c("genotypes", class(rez))
	if (!is.null(attr(a, "ped")))
		attr(rez, "ped") <- attr(a, "ped")
	if (!is.null(attr(a, "map")) & !is.null(attr(b, "map")))
		attr(rez, "map") <- rbind( attr(a, "map"), attr(b, "map") )
	
	return(rez)
	
}

## add a method to overload merge(), again keeping attributes
merge.genotypes <- function(a, b, join = c("inner","left"), ...) {
	
	if (!inherits(b, "genotypes"))
		stop("Only willing to merge two objects of class 'genotypes.'")
	
	if (length(intersect(colnames(a), colnames(b))))
		stop("Some samples are shared: that merge isn't implemented yet.")
	
	if (mode(a) != mode(b))
		stop("The two objects appear to have different modes (character vs. numeric).")
	
	if (!is.null(attr(a, "alleles"))) {
		if (!is.null(attr(b, "alleles"))) {
			if (attr(b, "alleles") != attr(a, "alleles"))
				warning("Alleles are coded differently in the two datasets; consider fixing that before merging.")
		}
	}
	
	o <- setNames( 1:nrow(a), rownames(a) )
	keep <- intersect(rownames(a), rownames(b))
	if (!length(keep))
		stop("No shared markers; result would be empty.")
	
	message(paste0("Set A has ", nrow(a), " markers x ", ncol(a), " samples."))
	message(paste0("Set B has ", nrow(b), " markers x ", ncol(b), " samples."))
	
	join <- match.arg(join)
	new.o <- keep[ order(o[keep]) ]
	if (join == "inner") {
		## keep intersection of marker sets
		rez <- cbind( unclass(a)[ new.o, ],
					  unclass(b)[ new.o, ] )
		if (.has.valid.intensity(a) && .has.valid.intensity(b)) {
			attr(rez, "intensity") <- list( x = cbind(attr(a, "intensity")$x[ new.o, ], attr(b, "intensity")$x[ new.o, ]),
											y = cbind(attr(a, "intensity")$y[ new.o, ], attr(b, "intensity")$y[ new.o, ]) )
			attr(rez, "normalized") <- .null.false(attr(a, "normalized")) && .null.false(attr(a, "normalized"))
		}
	}
	else
		stop("Not yet implemented: merges other than 'inner join'.")
	
	message(paste0("Merged set has ", nrow(rez), " markers x ", ncol(rez), " samples."))
	
	## add class info
	class(rez) <- c("genotypes", class(rez))
	
	## merge markers and family information
	if (!is.null(attr(a, "map")))
		attr(rez, "map") <- attr(a, "map")[ keep[ order(o[keep]) ], ]
	if (!is.null(attr(a, "ped")) & !is.null(attr(b, "ped")))
		attr(rez, "ped") <- rbind(attr(a, "ped"), attr(b, "ped"))
	
	return(rez)
	
}

## internal helpers for validating the 'genotypes' data structure and its parts

.is.valid.map <- function(map, ...) {
	return( all(colnames(map)[1:4] == c("chr","marker","cM","pos")) )
}

.has.valid.map <- function(gty, ...) {

	map <- attr(gty, "map")
	rez <- FALSE

	if (!is.null(map))
		rez <- .is.valid.map(map)
	if (is.na(rez))
		rez <- FALSE

	return(rez)

}

.has.valid.ped <- function(gty, ...) {
	
	## TODO
	return( !is.null(attr(gty, "ped")) )
	
}

## convert NULLs to FALSE (mostly in calls to attr(x,"not.here"))
.null.false <- function(x) {
	
	if (is.null(x))
		FALSE
	else x
	
}

.has.valid.intensity <- function(gty, ...) {
	
	rez <- FALSE
	if (!is.null(attr(gty, "intensity")))
		if (is.list(attr(gty, "intensity")) && length(attr(gty, "intensity")) == 2)
			if ( all(dim(attr(gty, "intensity")$x) == dim(attr(gty, "intensity")$y)) &&
				 	all(dim(attr(gty, "intensity")$x) == dim(gty)) )
				rez <- TRUE
	
	return(rez)
	
}

## apply a function over samples in a genotype matrix, by sample groups
gapply <- function(gty, expr, fn = NULL, unclass = FALSE, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Only willing to subset an object of class 'genotypes.'")
	
	e <- substitute(expr)
	if (!is.null(attr(gty, "ped")))
		r <- eval(e, attr(gty, "ped"), parent.frame())
	else
		r <- eval(e)
	
	if (unclass)
		gty <- unclass(gty)
	
	vals <- unique(r)
	rez <- lapply(vals, function(v) {
		fn(gty[ ,(r == v), drop = FALSE ], ...)
	})
	names(rez) <- vals
	return(rez)
	
}

maf <- function(x, ...) {
	rowMeans(x, na.rm = TRUE)/2
}

## convert calls from Hyuna Yang's mouseDivGeno to biallelic matrix
convert.mdg <- function(gty, map, vino = NA, ...) {

	if (is.null(rownames(map)))
		rownames(map) <- as.character(map[,2])

	message(paste0("Genotype matrix is ", nrow(gty), " markers x ", ncol(gty), " samples."))
	keep <- seq_len(nrow(gty))
	if (!is.null(rownames(gty))) {
		markers <- as.character(map[,2])
		keep <- intersect(markers, rownames(gty))
		map <- map[ keep, ]
		gty <- gty[ keep, ]
	}
	else {
		if (nrow(gty) != nrow(map))
			stop(paste0("Number of markers in map (", nrow(map), ") doesn't match dimension of genotype matrix."))
	}
	message(paste0("Keeping calls at ",length(keep), " markers."))
	vino.is <- c("A", "H", "B")[ vino ][1]
	message(paste0("VINO calls will be converted to '",vino.is,"'."))

	## initialize new matrix and set row/col names
	rez <- matrix(ncol = ncol(gty), nrow = length(keep))
	if (length(keep) == nrow(gty))
		rownames(rez) <- rownames(gty)
	else
		rownames(rez) <- keep
	colnames(rez) <- colnames(gty)

	## recode genotypes
	rez[ which(gty < 1 | gty > 4, arr.ind = TRUE) ] <- NA
	rez[ which(gty == 1, arr.ind = TRUE) ] <- 0
	rez[ which(gty == 2, arr.ind = TRUE) ] <- 1
	rez[ which(gty == 3, arr.ind = TRUE) ] <- 2
	rez[ which(gty == 4, arr.ind = TRUE) ] <- vino

	## sort by chromosome and position
	o <- with(map, order(chr, pos, cM))
	rez <- rez[ o, ]
	map <- map[ o, ]

	## add class information
	class(rez) <- c("genotypes", class(rez))
	attr(rez, "map") <- map

	## done
	return(rez)

}

#map.genotypes <- function(gty, ...) {
#	attr(gty, "map")
#}

#samples.genotypes <- function(gty, ...) {
#	attr(gty, "ped")
#}

.is.geno.matrix <- function(gty, ...)  {

	if (!(is.matrix(gty) & mode(gty) == "character")) {
		if (is.data.frame(gty)) {
			cols <- 3
			if ("cm" %in% tolower(colnames(gty)))
				cols <- 4
			message("Input was dataframe; assuming its first 3 (4) columns are chr-marker-(cM)-pos and converting remainder to character")
			rez <- as.matrix( gty[ ,-(1:cols) ] )
			attr(rez, "map") <- gty[ ,1:cols ]
			return(rez)
		}
		else
			stop("Please supply a character matrix (markers x samples)")
	}
	else
		return(gty)

}

filter.sites <- function(gty, maxn = 0.999999, maxhet = 0.999999, drop = TRUE, ...) {

	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")

	n <- ncol(gty)
	Ns <- apply(gty, 1, function(x) sum(x == "N" | is.na(x)))
	Hs <- apply(gty, 1, function(x) sum(x == "H" | x == 1))

	if (maxn < 1 && maxhet < 1) {
		Ns <- Ns/n
		Hs <- Hs/n
	}
	
	keep <- (Ns <= maxn & Hs <= maxhet)
	
	if (drop) {
		rez <- gty[ keep, ]
		attr(rez, "filtered") <- which(!keep)
	}
	else {
		rez <- keep
	}

	return(rez)

}

## fast way to find most common item in vector of numbers
top1 <- function(x, ...) {
	which.max( tabulate(x) )
}

## more friendly major-allele finder, using allele calls
major.allele <- function(x, nocall = c("N","H"), ...) {

	if (!is.character(x) & !is.factor(x)) {
		warning("Converting genotypes to character; this might not be what you wanted.")
		x <- as.character(x)
	}

	xx <- x[ !(x %in% nocall | is.na(x)) ]
	if (length(xx)) {
		return( names( sort(table(x), decreasing = TRUE) )[1] )
	}
	else {
		return(NA)
	}

}

minor.allele <- function(x, nocall = c("N","H"), ...) {

	if (!is.character(x) & !is.factor(x)) {
		warning("Converting genotypes to character; this might not be what you wanted.")
		x <- as.character(x)
	}

	maj <- major.allele(x, nocall = nocall)
	xx <- x[ !(x %in% nocall | is.na(x) | x == maj) ]
	if (length(xx)) {
		return( names( sort(table(x), decreasing = FALSE) )[1] )
	}
	else {
		return(NA)
	}

}

consensus.geno <- function(gty, ...) {

	warning("function consensus.geno() is deprecated!!")
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")

	if (ncol(gty) > 1)
		return( apply(gty, 1, major.allele, ...) )
	else
		## the trivial case: just return the input, as a vector
		return(as.vector(gty))

}

consensus <- function(gty, nas.allowed = 0.0, ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	rez <- apply(gty, 1, function(x) which.max(tabulate(x+1, nbins = 3))-1)
	nas <- apply(gty, 1, function(x) sum(is.na(x))/length(x) > nas.allowed)
	rez[nas] <- NA
	
	return(rez)
	#return(as.matrix(rez))
	
}

is.segregating <- function(x, ...) {

	if (is.character(x)) {
		x <- factor(x)
		if (any(!is.na(x)) | any(x != "N")) {
			flag <- sum(tabulate(x, nbins = 3) > 0) > 1
			return(flag | any(x == "H", na.rm = TRUE))
		}
		else
			return(FALSE)
	}
	else {
		if (any(!is.na(x))) {
			flag <- sum(tabulate(x+1, nbins = 3) > 0) > 1
			return(flag | any(x == 1, na.rm = TRUE))
		}
		else
			return(FALSE)
	}

}

segregating <- function(gty, ...) {

	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")

	keep <- apply(gty, 1, is.segregating)
	return(gty[ keep, ])

}

# is.fixed.diff <- function(gty, allowed = c("A","C","G","T"), ...) {
#
# 	if (!inherits(gty, "genotypes"))
# 		stop("Please supply an object of class 'genotypes.'")
#
# 	ns <- apply(gty, 1, function(x) any(is.na(x) | x == "N"))
# 	hs <- apply(gty, 1, function(x) any(x == 1 | x == "H"))
# 	diffs <- apply(gty, 1, function(x) {
# 		if (is.character(x))
# 			x <- as.numeric(factor(x, levels = allowed))
# 		xmin <- min(x, na.rm = TRUE)
# 		if (!is.finite(xmin))
# 			xmin <- 0
# 		if (min(x, na.rm = TRUE) < 1)
# 			x <- x + 1
# 		sum(tabulate(x) > 0) > 1
# 	})
#
# 	return(!ns & !hs & diffs)
#
# }

is.fixed.diff <- function(gty, ...) {

	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")

	if (!is.numeric(gty) || ncol(gty) != 2)
		stop("This function is for numeric genotypes for pairs of samples only.")

	gty <- unclass(gty)
	#print(head(gty))
	attr(gty, "map") <- NULL

	diffs <- (gty[ ,1 ] != gty[ ,2 ])
	ns <- is.na(rowSums(gty))
	hs <- (gty[ ,1 ] == 1) | (gty[ ,2 ] ==1)

	return(diffs & !ns & !hs)

}

heterozygosity <- function(x, het.char = "H", ...) {

	if (is.factor(x) || is.character(x))
		sum(x == het.char, na.rm = TRUE)/sum(!is.na(x))
	else
		sum(x == 1, na.rm = TRUE)/sum(!is.na(x))

}

recode.genotypes <- function(gty, mode = c("pass","01","native","relative"),
							 allowed = c("A","C","G","T","H"), alleles = NULL, ...) {

	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (is.null(attr(gty, "alleles"))) {
		attr(gty, "alleles") <- FALSE
	}
		
	## do the genotypes come with a map and reference alleles?
	alleles <- NULL
	if (.has.valid.map(gty))
		if (ncol(attr(gty, "map")) >= 6)
			alleles <- as.matrix(attr(gty, "map")[ ,5:6 ])

	## a suite of recoding functions which operate on (marker-wise) vectors of genotype calls

	# recode 0=major allele, 1=het, 2=minor allele
	.recode.numeric.by.freq <- function(calls, alleles = NULL) {

		calls <- factor( as.character(calls), levels = allowed )
		maj <- allowed[ top1(calls[ calls != "H" ]) ]
		new.calls <- rep(NA, length(calls))
		new.calls[ calls == "H" ] <- 1
		new.calls[ calls == maj ] <- 0
		new.calls[ is.na(new.calls) & !is.na(calls) ] <- 2
		return(new.calls)

	}

	# recode 0=A1, 1=het, 2=A2
	.recode.numeric.by.ref <- function(calls, alleles = NULL) {

		new.calls <- rep(NA, length(calls))
		new.calls[ calls == alleles[1] ] <- 0
		new.calls[ calls == alleles[2] ] <- 2
		new.calls[ calls == "H" ] <- 1
		return(new.calls)

	}

	# recode same as above, only backwards
	.recode.character.by.ref <- function(calls, alleles = NULL) {

		new.calls <- rep(NA, length(calls))
		new.calls[ calls == 0 ] <- as.character(alleles[1])
		new.calls[ calls == 2 ] <- as.character(alleles[2])
		new.calls[ calls == 1 ] <- "H"
		return(new.calls)

	}

	mode <- match.arg(mode)
	coding <- attr(gty, "alleles")
	recode.as <- FALSE
	if (mode == "pass")
		converter <- identity
	else if (mode == "01") {
		if ((!coding || coding == "01") && is.numeric(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			recode.as <- "01"
			converter <- identity
		}
		else if (!is.null(alleles)) {
			converter <- .recode.numeric.by.ref
			recode.as <- "01"
			message("Recoding to 0/1/2 using reference alleles.")
		}
		else {
			converter <- .recode.numeric.by.freq
			recode.as <- "relative"
			message("Recoding to 0/1/2 using empirical frequencies.")
		}
	}
	else if (mode == "relative") {
		if ((!coding || coding == "relative") && is.numeric(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			recode.as <- "relative"
			converter <- identity
		}
		else {
			converter <- .recode.numeric.by.freq
			recode.as <- "relative"
			message("Recoding to 0/1/2 using empirical frequencies.")
		}
	}
	else if (mode == "native")
		if (is.character(gty)) {
			message("Nothing to do; genotypes already in requested coding.")
			converter <- identity
		}
		else
			if (!is.null(alleles)) {
				converter <- .recode.character.by.ref
				recode.as <- "native"
			}
			else
				stop("Can only convert genotypes numeric->character given some reference alleles.")

	.gty <- unclass(gty)
	rez <- matrix(NA, nrow = nrow(.gty), ncol = ncol(.gty))
	colnames(rez) <- colnames(.gty)
	rownames(rez) <- rownames(.gty)
	for (i in seq_len(nrow(.gty))) {
		rez[ i, ] <- converter(.gty[i,], alleles = alleles[i,])
	}

	if (.has.valid.map(gty))
		attr(rez, "map") <- attr(gty, "map")
	if (.has.valid.ped(gty))
		attr(rez, "ped") <- attr(gty, "ped")
	
	for (a in c("intensity","normalized","filter.sites","filter.samples")) {
		if (!is.null(attr(gty, a)))
			attr(rez, a) <- attr(gty, a)
	}
	class(rez) <- c("genotypes", class(rez))
	attr(rez, "alleles") <- recode.as
	return(rez)

}

## code genotypes as sharing 0,1,2 alleles IBS with a reference sample.
ibs.to.reference <- function(gty, ref = 1, ...) {

	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")

	if (!is.numeric(gty))
		gty <- recode.genotypes(gty, mode = "01")

	if (length(ref) == 1) {
		## assume ref points to the reference column
		ref <- unclass(gty)[ ,ref, drop = TRUE ]
	}
	else if (length(ref) == nrow(gty) && mode(gty) == mode(ref)) {
		## assume ref is a vector of reference genotypes
		ref <- ref
	}
	else {
		stop("Can't understand how reference sample (or genotypes) was specified.")
	}

	.ibs.to <- function(x, y) {
		rez <- numeric(length(x))
		rez <- 2-abs(x-y)
		rez[ y == 1 ] <- 1
		return(rez)
	}

	rez <- apply(unclass(gty), 2, .ibs.to, ref)
	attributes(rez) <- attributes(gty)
	return(rez)

}

## convert genotypes (in matrix or datafame) to plink's tped format (markers x samples)
## <gty> = genotype matrix (markers x samples), or dataframe with map (see below) + genotypes
## <map> = dataframe with up to 6 columns: chr, marker, cM pos, bp pos, allele 1, allele 2
geno.to.tped <- function(gty, map = NULL, alleles = NULL, nocalls = c("N"), het = "H", resort = FALSE, ...) {

	.expand.numeric.geno <- function(col, ...) {
		rez <- matrix(0, ncol = 2, nrow = length(col))
		rez[ col == 0, 1 ] <- rez[ col == 0, 2 ] <- 1
		rez[ col == 2, 1 ] <- rez[ col == 2, 2 ] <- 2
		rez[ col == 1, 1 ] <- 1
		rez[ col == 1, 2 ] <- 2
		return(rez)
	}

	.expand.character.geno <- function(col, ...) {
		rez <- matrix("0", ncol = 2, nrow = length(col))
		for (a in 1:2) {
			rez[ col == alleles[,a],1 ] <- rez[ col == alleles[,a],2 ] <- alleles[ col == alleles[,a],a ]
		}
		rez[ col == "H",1 ] <- alleles[ col == "H",1 ]
		rez[ col == "H",2 ] <- alleles[ col == "H",1 ]
		#print(dim(rez))
		return(rez)
	}

	if (is.matrix(gty)) {
		message("Assuming genotypes supplied as matrix.")
		if (is.null(map)) {
			if ("map" %in% names(attributes(gty)))
				map <- attr(gty, "map")
			else
				stop("If supplying genotypes as matrix, must also supply a map.")
		}
	}
	else if (is.data.frame(gty)) {
		message("Assuming genotypes supplied as dataframe including map.")
		if (is.null(map))
			map <- gty[ ,1:4 ]
		gty <- as.matrix(gty[ ,-(1:4) ])
	}
	else
		stop("Must supply genotypes as either dataframe (including map) or matrix+map.")

	if (nrow(gty) != nrow(map))
		stop("Dimensions of genotype matrix and map don't match.")

	## set up genotype conversion
	converter <- function() {}
	if (is.numeric(gty))
		converter <- .expand.numeric.geno
	else {
		converter <- .expand.character.geno
		if (is.null(alleles)) {
			if (ncol(map) < 6)
				stop("Need to specify alleles for character genotypes; might not discover them from the data.")
			else {
				alleles <- as.matrix(map[ ,5:6 ])
				rownames(alleles) <- as.character(map[,2])
			}
		}
		else {
			alleles <- alleles[ as.character(map$marker), ]
		}
	}

	## convert chromosome names, if needed
	if (any(grepl("^chr", map[,1]))) {
		message("Converting chromosome names to plink style.")
		map[,1] <- convert.names(map[,1], to = "plink")
	}

	## convert 1-column genoypes to plink's 2-column format (TODO: what about phasing?)
	new.geno <- matrix(ncol = 2*ncol(gty), nrow = nrow(gty))
	pb <- txtProgressBar(min = 0, max = ncol(gty), style = 3)
	for (i in 1:ncol(gty)) {
		new.geno[ ,(2*i-1):(2*i) ] <- converter(gty[,i])
		setTxtProgressBar(pb, i)
	}
	## add marker map
	rez <- data.frame(map[,1:4], new.geno)
	if (resort)
		rez <- rez[ order(rez[,1], rez[,3], rez[,4]), ]

	## zero out unknown positions (genetic or physical)
	rez[ is.na(rez[,3]),3 ] <- 0
	rez[ is.na(rez[,4]),4 ] <- 0

	## strip markers which plink will not like
	discard <- !(rez[,1] %in% c(1:19,"X","Y","M","MT"))
	if (sum(discard) > 0) {
		rez <- rez[ !discard, ]
		warning(paste("Discarded", sum(discard), "markers on unknown chromosomes;", sum(!discard), "markers remain."))
	}

	## bless object and return
	class(rez) <- c("tped", class(rez))
	return(rez)

}

## replace an existing map with a new one, possibly with allele info
## if <newmap> has rownames, they should map old marker names to new ones
replace.map <- function(gty, newmap, ...) {

	if (!inherits(gty, "genotypes") | !.has.valid.map(gty))
		stop("Please supply a genotypes object with map as an attribute.")

	map <- attr(gty, "map")
	map$order <- 1:nrow(map)
	rownames(map) <- as.character(map$marker)
	message(paste("Starting with", nrow(map), "markers..."))

	if (!is.null(rownames(newmap))) {
		## extract only markers which are (1) present in current map; (2) accounted for in new map
		m <- match(as.character(map$marker), rownames(newmap), nomatch = 0)
		gty <- gty[ which(m > 0), ]
		rownames(gty) <- newmap[ m[m > 0],"marker" ]
		attr(gty, "map") <- newmap[ m[m > 0], ]
	}
	else {
		attr(gty, "map") <- newmap
	}

	## sort the result
	o <- order( attr(gty, "map")$chr, attr(gty, "map")$pos, attr(gty, "map")$marker )
	map <- attr(gty, "map")
	gty <- gty[ o, ]
	
	## sort intensities too, if present
	if (.has.valid.intensity(gty)) {
		if (is.list(attr(gty, "intensity")) && length(attr(gty, "intensity")) == 2) {
			x.new <- attr(gty, "intensity")$x[ o, ]
			y.new <- attr(gty, "intensity")$y[ o, ]
			attr(gty, "intensity") <- list(x = x.new, y = y.new)
		}
	}
	
	cols <- colnames(map)
	required <- c("chr","marker","cM","pos")
	attr(gty, "map") <- map[ o,c(required, setdiff(cols, required)) ]
	
	message(paste("...and ending with", nrow(attr(gty, "map")), "markers."))
	if (nrow(attr(gty, "map")) != nrow(gty))
		stop("Map no longer matches genotypes object.")

	return(gty)

}

## utility function to convert map+genotype matrix (columns marker-chr-pos-{genotypes}) to R/qtl's "CSVs" input format
## gty (or map, if supplied separately) should have rownames which are marker IDs
as.rqtl.geno <- function(gty, map = NULL, ...) {

	if (is.matrix(gty)) {
		if (is.null(map))
			stop("If supplying genotypes as matrix, must also supply a map.")
	}
	else if (is.data.frame(gty)) {
		gty <- as.matrix(gty[ ,-(1:3) ])
		if (is.null(map))
			map <- gty[ ,1:3 ]
	}
	else
		stop("Must supply genotypes as either dataframe (including map) or matrix+map.")
	if (is.null(rownames(map)))
	rownames(map) <- as.character(map[,1])

	rez <- matrix("", ncol = nrow(gty) + 1, nrow = ncol(gty) + 2)
	colnames(rez) <- c("id", rownames(map))

	rez[ 1,2:ncol(rez) ] <- map[,2]
	rez[ 2,2:ncol(rez) ] <- map[,3]
	rez[ 3:nrow(rez),2:ncol(rez) ] <- t(gty)
	rez[ 3:nrow(rez),1 ] <- colnames(gty)

	return(rez)

}

## read tped(-like) file and make it a dataframe
## tped format: <chr> <marker> <cM> <pos> (<allele 1> <allele 2>)n
## biallelic markers assumed
read.tped <- function(..., samples = NULL) {

	df <- read.table(..., stringsAsFactors = FALSE)
	colnames(df)[1:4] <- c("chr","marker","cM","pos")
	newdf <- df[,1:4]
	as.char <- FALSE
	if (is.character(df[,5]) | is.factor(df[,5]))
		as.char <- TRUE

	for (i in 0:((ncol(df)-4)/2-1)) {
		is.het <- df[ ,5+2*i ] != df[ ,5+2*i+1 ]
		if (as.char) {
			newcol <- as.character(df[ ,5+2*i ])
			newcol[is.het] <- "H"
		}
		else {
			newcol <- df[ ,5+2*i ]
			newcol[ newcol == 0 ] <- NA
			newcol <- (newcol-1)*2
			newcol[is.het] <- 1
		}
		newdf <- cbind(newdf, newcol)
	}

	if (length(samples) == ncol(newdf)-4)
		colnames(newdf)[ 5:ncol(newdf) ] <- samples
	return(newdf)

}

.chr.greedy.ibd <- function(gty, mismatches = 0, ...) {

	gty <- unclass(gty)

	segstart <- NA
	segend <- NA
	n <- 0
	mm <- 0
	segments <- matrix(nrow = 0, ncol = 4)
	for (i in 1:nrow(gty)) {

		## ignore no-calls
		if (any(gty[ i, ] == "N" | is.na(gty[ i, ])))
			next

		## het sites are ambiguous; be optimistic and assume the segment extends
		if (any(gty[ i, ] == "H" | gty[ i, ] == 1)) {
			if (is.na(segstart)) {
				segstart <- segend <- i
			}
			else
				segend <- i
			n <- n+1
		}
		## IBS: segment definitely extends
		else if (gty[ i,1 ] == gty[ i,2 ]) {
			if (is.na(segstart)) {
				segstart <- segend <- i
			}
			else
				segend <- i
			n <- n+1
		}
		## definitely not IBS: terminate segment
		else if (gty[ i,1 ] != gty[ i,2 ]) {
			if (!is.na(segstart)) {
				if (mm == mismatches) {
					if (n > 1) {
						segments <- rbind(segments, c(segstart, segend, n, 1-(mm/n)))
					}
					segstart <- NA
					segend <- NA
					n <- 0
					mm <- 0
				}
				else {
					segend <- i
					n <- n+1
					mm <- mm+1
				}
			}
		}

	} #end for

	if (!is.na(segstart)) {
		segend <- nrow(gty)
		if (n > 1) {
			segments <- rbind(segments, c(segstart, segend, n, 1-(mm/n)))
		}
	}

	return(segments)

}

## use greedy algorithm to define maximal intervals of pairwise identity-by-descent
## <gty> = a matrix (character or numeric) of genotypes, assumed sorted by chr and pos
## <map> = dataframe of marker,chr,pos,[cM] specifying marker locations
## <pairs> = 2-col matrix of pairs over which to compute IBD, one row per pair
greedy.ibd <- function(gty, map = NULL, pairs = t(combn(ncol(gty), 2)), mismatches = 0, ...) {

	require(plyr)

	if (!inherits(gty, "genotypes"))
		gty <- .is.geno.matrix(gty)

	if (is.null(map))
		if (.has.valid.map(gty))
			map <- attr(gty, "map")
		else
			stop("Please supply a map, either as attribute of <gty> or as argument to this call.")

	segments <- data.frame(chr = character(), start = numeric(), end = numeric(),
						   p1 = character(), p2 = character(), ident = numeric())
	map$chr <- factor(map$chr)
	for (c in levels(map$chr)) {
		i <- (map$chr == c & !is.na(map$pos))
		cat("--- chromosome:", c, "(", sum(i, na.rm = TRUE), "markers )---\n")
		for (p in 1:nrow(pairs)) {
			cat("\t pair (", pairs[p,1], ",", pairs[p,2], ")\n")
			rez <- .chr.greedy.ibd(gty[i, pairs[p,]], mismatches)
			if (nrow(rez) > 0) {
				#if (any(is.na(rez[,2])))
				#	rez[ is.na(rez[,2]),2 ] <- max(map$pos[i], na.rm = TRUE)
				segments <- rbind(segments,
								  data.frame(chr = c, start = map$pos[ which(i)[rez[,1]] ],
								  		   end = map$pos[ which(i)[rez[,2]] ],
								  		   n = rez[,3], p1 = pairs[p,1], p2 = pairs[p,2], ident = rez[,4]))
			}
		}
	}

	if (!is.character(pairs) & !is.null(colnames(gty))) {
		segments$p1 <- colnames(gty)[ segments$p1 ]
		segments$p2 <- colnames(gty)[ segments$p2 ]
	}

	return(segments)

}

# ibd.hmm <- function(...) {
#
# }

## find IBD segments using an HMM, per chromosome
## HMM states: IBS0, IBS1, IBS2
## HMM symbols = # alleles shared: 0, 1, 2
## transitions are functions of map distances
# .chr.ibd.hmm <- function(geno, map = NULL, eps = 0.6, Ne = 100, ...) {
#
# 	## check that genotypes and map are present; strip markers with missing metadata
# 	if (is.null(map))
# 		if (!is.null(attr(geno, "map")))
# 			map <- attr(geno, "map")
# 	nas <- attr(na.omit(map), "na.action")
# 	if (!is.null(nas) | length(nas)) {
# 		map <- map[ -nas, ]
# 		geno <- geno[ -nas, ]
# 	}
#
# 	## calculate per-interval transition probabilities using genetic map
# 	map$dist <- c(0, diff(map$cM))
# 	map$dist[ c(1, which(map$chr[ 2:(nrow(map)) ] != map$chr[ 1:(nrow(map)-1) ])+1) ] <- 0
# 	tmat <- array(0, dim = c(nrow(map), 3, 3))
# 	for (i in seq_len(nrow(map))) {
# 		this.t <- diag(3)
# 		d <- map$dist[i]/100
# 		r <- (4*Ne*d + 1)^(-1)
# 		this.t <- matrix(c( 1-r-r^2, r, r^2,
# 							r, 1-2*r, r,
# 							r^2, r, 1-r-r^2),
# 						 byrow = TRUE, nrow = 3)
# 		tmat[i,,] <- this.t
# 	}
#
# 	print(tmat[2,,])
#
# 	## emission probabilities: should scale this for MAF or something...
# 	emit <- eps*diag(3)
# 	emit[ upper.tri(emit)|lower.tri(emit) ] <- (1-eps)/2
#
# 	## prepare the hmm object
# 	hmm <- init.hmm(nrow(map), tmat, emit)
#
# 	## observed states: # alleles shared IBD
# 	## NB: heterozygous cases, as written, may only work if one sample is haploid/inbred
# 	sharing <- 1 # default: 1 allele shared IBD
# 	sharing[ geno[,1] == geno[,2] & geno[,1] != 1 ] <- 2 # could be a problem in het case
# 	sharing[ geno[,1] != geno[,2] & (geno[,1] == 1 | geno[,2] == 1) ] <- 1
# 	sharing[ geno[,1] != geno[,2] & !(geno[,1] == 1 | geno[,2] == 1) ] <- 0
# 	sharing[ is.na(sharing) ] <- 1
# 	#print(table(sharing, useNA = "always"))
# 	rez <- viterbi(hmm, sharing)
# 	return(rez$decoded)
#
# }

## convert bp->cM using linear interpolation on a reference map
## <x> = dataframe-like object with (at least) 2 columns: chr, pos
## <map> = plink-style genetic map as dataframe-like object (4 columns: chr, marker, cM, pos)
interpolate.map <- function(x, map, min.bp = 3e6, what = c("cM","bp"), seqlengths = NULL, ...) {

	require(plyr)

	x <- as.data.frame(x)
	if (ncol(x) < 2)
		stop("Need at least 2 columns of data (chr, pos)")

	if (!all(c("chr","pos") %in% colnames(x))) {
		warning("Assuming first 2 columns are (chr,pos)")
		colnames(x)[1:2] <- c("chr","pos")
	}
	x <- subset(x, pos > 0)
	x <- x[ with(x, order(chr, pos)), ]

	what <- match.arg(what)
	if (what == "cM")
		map <- map[ with(map, order(chr, pos, cM)), ]
	else
		map <- map[ with(map, order(chr, cM, pos)), ]
	x$chr <- factor(x$chr, levels = levels(map$chr))
	#print(head(x))
	
	.bp.to.cM <- function(chr, p, ref, z = min.bp, ...) {

		ref <- ref[ ref[,2] > 0, ]
		if (min(ref[,1], na.rm = TRUE) > min.bp)
			ref <- rbind(c(min.bp, 0), ref)
		i <- pmin(pmax(findInterval(p, ref[,1]), 1, na.rm = TRUE), nrow(ref)-1, na.rm = TRUE)
		bp <- p - ref[i,1]
		w <- (ref[i+1,1] - ref[i,1])
		offset <- ref[i,2]
		if (!is.null(seqlengths)) {
			maxlen <- seqlengths[ as.character(chr) ]
			overs <- p > max(ref[,1])
			bp[overs] <- p[overs] - max(ref[,1])
			offset[overs] <- max(ref[,2])
			w[overs] <- maxlen - bp[overs]
		}
		
		cm <- offset + (bp/w)*(ref[i+1,2] - ref[i,2])
		return(cm)

	}

	rez <- ddply(x, .(chr), function(cc) {
		#print(summary(cc))
		.map <- subset(map, chr == cc$chr[1] & !is.na(pos) & !is.na(cM))
		if (nrow(.map))
			if (what == "cM")
				data.frame(cc, cM = .bp.to.cM(cc$chr[1], cc$pos, .map[ ,c("pos","cM") ]))
			else
				data.frame(cc, bp = .bp.to.cM(cc$chr[1], cc$cM, .map[ ,c("cM","pos") ]))
		else
			if (what == "cM")
				data.frame(cc, cM = NA)
			else
				data.frame(cc, pos = NA)
	}, .progress = "none")

	rez <- rez[ with(rez, order(chr, pos, cM)), ]
	return(rez)

}

is.autosome <- function(x, ...) {

	as.character(x) %in% c(1:19, paste0("chr",1:19))

}

## return just autosomes from a genotype matrix or dataframe
autosomes <- function(gty, ...) {

	chr.col <- c("chr","chrm","chrom","seqnames")
	if (inherits(gty, "genotypes") & .has.valid.map(gty)) {
		i <- is.autosome(attr(gty, "map")[,1])
		return( gty[i,] )
	}
	else if (is.data.frame(gty) & any(chr.col %in% colnames(gty))) {
		col <- which(colnames(gty) %in% chr.col)[1]
		i <- is.autosome(gty[,col])
		return( gty[i,] )
	}
	else if (inherits(gty, "GRanges")) {
		return( gty[ is.autosome(as.vector(seqnames(gty))) ] )
	}
	else if (inherits(gty, "Seqinfo")) {
		chroms <- c(1:19)
		if (any(grepl("^chr", seqnames(gty))))
			chroms <- paste0("chr", chroms)
		return( gty[ chroms ] )
	}
	else {
		stop("Please supply genotype matrix with attached map, or a dataframe with column 'chr'.")
	}

}

normalchrs <- function(gty, ...) {
	
	.is.normal.chr <- function(y) as.character(y) %in% c(1:19, paste0("chr",1:19), "X", "chrX")
	
	chr.col <- c("chr","chrm","chrom","seqnames")
	if (inherits(gty, "genotypes") & .has.valid.map(gty)) {
		i <- .is.normal.chr(attr(gty, "map")[,1])
		return( gty[i,] )
	}
	else if (is.data.frame(gty) & any(chr.col %in% colnames(gty))) {
		col <- which(colnames(gty) %in% chr.col)[1]
		i <- .is.normal.chr(gty[,col])
		return( gty[i,] )
	}
	else if (inherits(gty, "GRanges")) {
		return( gty[ .is.normal.chr(as.vector(seqnames(gty))) ] )
	}
	else if (inherits(gty, "Seqinfo")) {
		chroms <- c(1:20)
		if (any(grepl("^chr", seqnames(gty))))
			chroms <- paste0("chr", chroms)
		return( gty[ chroms ] )
	}
	else {
		stop("Please supply genotype matrix with attached map, or a dataframe with column 'chr'.")
	}
	
}

## make a fictional mutliple sequence alignment from genotype calls
## 	by converting Ns to random non-allelic nucleotides
spoof.msa <- function(gty, mode = c("nt","character"), ...) {
	
	require(ape)
	
	if (!inherits(gty, "genotypes") | !is.character(gty) | !all(c("A1","A2") %in% colnames(attr(gty, "map"))))
		stop("Please supply 'genotypes' object, coded as characters, with map with allele information.")
	
	mode <- match.arg(mode)
	if (mode == "nt") {
		
		nts <- c("A","C","G","T")
		recoded <- matrix("H", nrow = nrow(gty), ncol = ncol(gty))
		for (i in seq_len(nrow(gty))) {
			alleles <- toupper( as.vector(unlist( attr(gty, "map")[ i,c("A1","A2") ] )) )
			alts <- sample(setdiff(nts, alleles), 2)
			row <- as.vector(gty[ i, ])
			row[ row == "N" ] <- alts[1]
			row[ row == "H" ] <- alts[2]
			recoded[ i, ] <- row
		}
		recoded[ recoded == "H" ] <- "N"
		
	}
	else if (mode == "character") {
		
		nts <- c("A","C","G","T")
		recoded <- matrix(0, nrow = nrow(gty), ncol = ncol(gty))
		for (i in seq_len(nrow(gty))) {
			alleles <- toupper( as.vector(unlist( attr(gty, "map")[ i,c("A1","A2") ] )) )
			alts <- c("H","N")
			old <- c(alleles, alts)
			vals <- seq_along(old)-1
			row <- as.vector(gty[ i, ])
			newrow <- numeric(length(row))
			for (j in seq_along(old)) {
				newrow[ row == old[j] ] <- vals[j]
			}
			recoded[ i, ] <- newrow
		}
		
	}
	
	rez <- as.alignment( t(recoded) )
	rez$nam <- colnames(gty)
	return(rez)
	
}