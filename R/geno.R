## utility functions for handling genotype matrices

## convert between styles of chromosome names (I prefer UCSC)
convert.names <- function(chrs, to = c("ncbi","ensembl","ucsc","plink"), muga.style = FALSE, ...) {
	
	ucsc <- paste0("chr", c(1:19, "X","Y","M"))
	ncbi <- c(1:19,"X","Y","MT")
	
	converter <- function()
		to <- match.arg(to)
	if (to == "ncbi" | to == "ensembl") {
		converter <- function(x) setNames(ncbi, ucsc)[x]
	}
	else if (to == "ucsc") {
		converter <- function(x) setNames(ucsc, ncbi)[x]
	}
	else if (to == "plink") {
		converter <- function(x) {
			new <- gsub("^chr","", x)
			new <- gsub("^M$","MT", x)
			return(new)
		}
	}
	
	new <- converter( as.character(chrs) )
	if (muga.style)
		new <- gsub("^MT$","M", chrs)
	
	return(new)
	
}

geno.to.matrix <- function(gty, ...) {
	
	require(reshape2)
	
	gty.mat <- dcast(gty, marker + chr + pos ~ sid, value.var = "call")
	gty.mat <- gty.mat[ with(gty.mat, order(chr, pos)), ]
	rownames(gty.mat) <- gty.mat$marker
	return(gty.mat)
	
}

.is.geno.matrix <- function(gty, ...)  {
	
	if (!(is.matrix(gty) & mode(gty) == "character")) {
		if (is.data.frame(gty)) {
			warning("Input was dataframe; assuming its first 3 columns are marker-chr-pos and converting remainder to character")
			return( as.matrix( gty[ ,-(1:3) ] ) )
		}
		else
			stop("Please supply a character matrix (markers x samples)")
	}
	else
		return(gty)

}

filter.sites <- function(gty, maxn = 1, maxhet = 1, ...) {
	
	gty <- .is.geno.matrix(gty)	
	
	n <- ncol(gty)
	Ns <- apply(gty, 1, function(x) sum(x == "N"))
	Hs <- apply(gty, 1, function(x) sum(x == "H"))
	
	keep <- (Ns/n <= maxn & Hs/n <= maxhet)
	rez <- gty[ keep, ]
	attr(rez, "filtered") <- which(!keep)
	
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
	
	gty <- .is.geno.matrix(gty)
	if (ncol(gty) > 1)
		return( apply(gty, 1, major.allele, ...) )
	else
		## the trivial case: just return the input, as a vector
		return(as.vector(gty))
	
}

recode.genotypes <- function(gty, alleles = c("A","C","G","T","H"), ...) {
	
	gty <- .is.geno.matrix(gty)
	
	.recode.site <- function(calls) {
		
		calls <- factor( as.character(calls), levels = alleles )
		maj <- alleles[ top1(calls[ calls != "H" ]) ]
		new.calls <- rep(NA, length(calls))
		new.calls[ calls == "H" ] <- 1
		new.calls[ calls == maj ] <- 0
		new.calls[ is.na(new.calls) & !is.na(calls) ] <- 2
		return(new.calls)
		
	}
	
	rez <- t( apply(gty, 1, .recode.site) )
	colnames(rez) <- colnames(gty)
	rownames(rez) <- rownames(gty)
	return(rez)
	
}

## convert genotypes (in matrix or datafame) to plink's tped format (markers x samples)
## <gty> = genotype matrix (markers x samples), or dataframe with map (see below) + genotypes
## <map> = dataframe with 4 columns: chr, marker, cM pos, bp pos
geno.to.tped <- function(gty, map = NULL, nocalls = c("N"), het = "H", ...) {
	
	.expand.numeric.geno <- function(col, ...) {
		rez <- matrix(0, ncol = 2, nrow = length(col))
		rez[ col == 0, 1 ] <- rez[ col == 0, 2 ] <- 1
		rez[ col == 2, 1 ] <- rez[ col == 2, 2 ] <- 2
		rez[ col == 1, 1 ] <- 1
		rez[ col == 1, 2 ] <- 2
		return(rez)
	}
	
	.expand.character.geno <- function(col, ...) {
		rez <- matrix("N", ncol = 2, nrow = length(col))
		called <- any(!(col %in% nocalls))
		alleles <- unique( col[ !(col %in% c(nocalls, het)) ] )
		if (length(alleles) & called) {
			## only one homozygote class seen; spoof other allele
			if (length(alleles) == 1)
				alleles <- c(alleles, "Z")
			for (a in alleles) {
				rez[ col == a,1 ] <- rez[ col == a,2 ] <- a
			}
			rez[ col == "H",1 ] <- alleles[1]
			rez[ col == "H",2 ] <- alleles[2]
		}
		print(dim(rez))
		return(rez)
	}
	
	if (is.matrix(gty)) {
		if (is.null(map))
			stop("If supplying genotypes as matrix, must also supply a map.")
	}
	else if (is.data.frame(gty)) {
		gty <- as.matrix(gty[ ,-(1:4) ])
		if (is.null(map))
			map <- gty[ ,1:4 ]
	}
	else
		stop("Must supply genotypes as either dataframe (including map) or matrix+map.")
	
	## set up genotype conversion
	converter <- function() {}
	if (is.numeric(gty))
		converter <- .expand.numeric.geno
	else
		converter <- .expand.character.geno
	
	## convert chromosome names, if needed
	if (any(grepl("^chr", map[,1])))
		map[1,] <- convert.names(map[1,], to = "plink")
	
	## convert 1-column genoypes to plink's 2-column format (TODO: what about phasing?)
	new.geno <- matrix(ncol = 0, nrow = nrow(gty))
	for (i in 1:ncol(gty)) {
		new.geno <- cbind(new.geno, converter(gty[,i]))
	}
	## add marker map
	rez <- data.frame(map, new.geno)
	
	return(rez)
	
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
	
	df <- read.table(...)
	colnames(df)[1:4] <- c("chr","marker","cM","pos")
	newdf <- df[,1:4]
	for (i in 0:((ncol(df)-4)/2-1)) {
		is.het <- df[ ,5+2*i ] != df[ ,5+2*i+1 ]
		newcol <- as.character(df[ ,5+2*i ])
		newcol[is.het] <- "H"
		newdf <- cbind(newdf, newcol)
		
	}
	
	if (length(samples) == ncol(newdf)-4)
		colnames(newdf)[ 5:ncol(newdf) ] <- samples
	return(newdf)
	
}