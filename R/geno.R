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

geno.to.matrix <- function(gty, map = NULL, ...) {
	
	require(reshape2)

	## if map specified, attach it to get cM positions
	if (!is.null(map)) {
		nsnps.before <- length(unique(gty$marker))
		message(paste("Attaching map: started with", nsnps.before, "markers"))
		gty <- merge(gty, map[ ,c("marker","cM") ])
		nsnps.after <- length(unique(gty$marker))
		message(paste("Done attaching map: ended with", nsnps.after, "markers"))
	}
	
	## if genetic position included, keep it
	fm <- "chr + marker + pos ~ sid"
	if ("cM" %in% colnames(gty))
		fm <- "chr + marker + cM + pos ~ sid"
	
	gty.mat <- dcast(gty, as.formula(fm), value.var = "call")
	gty.mat <- gty.mat[ with(gty.mat, order(chr, pos)), ]
	rownames(gty.mat) <- gty.mat$marker
	return(gty.mat)
	
}

.is.geno.matrix <- function(gty, ...)  {
	
	if (!(is.matrix(gty) & mode(gty) == "character")) {
		if (is.data.frame(gty)) {
			cols <- 3
			if ("cm" %in% tolower(colnames(gty)))
				cols <- 4
			message("Input was dataframe; assuming its first 3 (4) columns are chr-marker-(cM)-pos and converting remainder to character")
			return( as.matrix( gty[ ,-(1:4) ] ) )
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
		if (is.null(map))
			map <- gty[ ,1:4 ]
		gty <- as.matrix(gty[ ,-(1:4) ])
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
	class(rez) <- c("tped", class(rez))
	
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

.chr.greedy.ibd <- function(gty, pair, ...) {
	
	segstart <- NA
	segend <- NA
	n <- 0
	segments <- matrix(nrow = 0, ncol = 3)
	for (i in 1:nrow(gty)) {
		
		## ignore no-calls
		if (any(gty[i,pair] == "N"))
			next
		
		## het sites are ambiguous; be optimistic and assume the segment extends
		if (any(gty[ i,pair ] == "H")) {
			if (is.na(segstart))
				segstart <- segend <- i
			else
				segend <- i
			n <- n+1
		}
		
		## IBS: segment definitely extends
		if (gty[ i,pair[1] ] == gty[ i,pair[2] ]) {
			if (is.na(segstart))
				segstart <- segend <- i
			else
				segend <- i
			n <- n+1
		}
		## definitely not IBS: terminate segment
		else {
			if (!is.na(segstart)) {
				segments <- rbind(segments, c(segstart, segend, n))
				segstart <- NA
				segend <- NA
				n <- 0
			}
		}
		
	} #end for
	
	return(segments)
	
}

## use greedy algorithm to define maximal intervals of pairwise identity-by-descent 
## <gty> = a matrix (character or numeric) of genotypes, assumed sorted by chr and pos
## <map> = dataframe of marker,chr,pos,[cM] specifying marker locations
## <pairs> = 2-col matrix of pairs over which to compute IBD, one row per pair
greedy.ibd <- function(gty, map, pairs = t(combn(ncol(gty), 2)), ...) {
	
	require(plyr)
	
	if (!is.matrix(gty))
		gty <- as.matrix(gty)
	
	segments <- data.frame(chr = character(), start = numeric(), end = numeric(), p1 = numeric(), p2 = numeric())
	map$chr <- factor(map$chr)
	for (c in levels(map$chr)) {
		cat("--- chromosome:", c, "---\n")
		i <- (map$chr == c)
		for (p in 1:nrow(pairs)) {
			cat("\t pair (", pairs[p,1], ",", pairs[p,2], ")\n")
			rez <- .chr.greedy.ibd(gty[i,], pairs[p,])
			print(head(rez))
			segments <- rbind(segments,
							  data.frame(chr = c,
							  		   start = map$pos[ which(i)[rez[,1]] ], end = map$pos[ which(i)[rez[,2]] ],
							  		   n = rez[,3], p1 = pairs[p,1], p2 = pairs[p,2]))
		} 
	}
	
	if (!is.null(colnames(gty))) {
		segments$p1 <- colnames(gty)[ segments$p1 ]
		segments$p2 <- colnames(gty)[ segments$p2 ]
	}
	
	return(segments)
	
}

## create a fake tfam-formatted pedigree file for plink
spoof.tfam <- function(ids, sex = 0, ...) {
	
	resex <- 0
	resex[ grepl("^[mM]", sex) ] <- 1
	resex[ grepl("^[fF]", sex) ] <- 2
	
	fam <- data.frame(fid = seq_along(ids), id = seq_along(ids),
					  pid = length(ids)+1, mid = length(ids)+2, sex = resex, phe = -9)
	return(fam)
	
}

## shortcut to write tab-separated output files without rownames, colnames, quotes
write.plink.file <- function(...) {
	write.table(..., col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

## compute genome-wide IBD estimate (\hat{pi}) with plink, optionally with initial LD pruning step
ibd.plink <- function(gty, map = NULL, ped = NULL, where = tempdir(), prefix = "stuff",
					  flags = "--nonfounders", prune = FALSE, prune.by = "--make-founders --indep 50 5 2", ...) {
	
	require(Matrix)
	
	## check if plink installed
	rez <- system("plink --version")
	if (rez != 0) {
		stop("Can't find plink executable; is it in your $PATH?")
	}
	
	## convert input to tped, if it isn't
	if (!inherits(gty, "tped"))
		gty <- geno.to.tped(gty, map)
	
	nind <- sum(!(tolower(colnames(gty)) %in% c("chr","pos","cm","id","marker")))/2
	gty <- subset(gty, !is.na(cM) & !is.na(pos) & pos > 0)
	
	if (is.null(ped))
		ped <- spoof.tfam(1:nind)
	
	prefix <- file.path(where, prefix)
	system(paste0("rm ", prefix, "*"))
	
	write.plink.file(ped, paste0(prefix, ".tfam"))
	write.plink.file(gty, paste0(prefix, ".tped"))
	
	cmd <- paste("plink --tfile", prefix, "--make-bed --out", prefix)
	system(cmd, intern = FALSE)
	
	if (prune) {
		prefix.new <- paste0(prefix, ".pruned")
		cmd <- paste("plink --bfile", prefix, prune.by, "--out", prefix.new)
		system(cmd, intern = FALSE)
		cmd <- paste("plink --bfile", prefix, "--extract", paste0(prefix.new, ".prune.in"), "--make-bed --out", prefix.new)
		system(cmd, intern = FALSE)
		if ( !file.exists(paste0(prefix.new, ".prune.in")) ) {
			stop( paste("LD pruning failed; command was '", cmd,"'") )
		}
		else {
	
		}
		prefix <- prefix.new
	}
	
	cmd <- paste("plink --bfile", prefix, "--genome --out", prefix.new, flags)
	rez <- system(cmd, intern = FALSE)
	ibd.file <- paste0(prefix.new, ".genome")
	if (file.exists(ibd.file)) {
		ibd <- read.table(ibd.file, header = TRUE, strip.white = TRUE)
		K <- matrix(NA, nrow = nind, ncol = nind)
		for (i in 1:nrow(ibd)) {
			K[ ibd$IID1[i], ibd$IID2[i] ] <- ibd$PI_HAT[i]
		}
		return(Matrix(K, sparse = TRUE))
	}
	else {
		stop( paste0("Call to plink failed; command was: '",cmd,"'") )
	}
	
}