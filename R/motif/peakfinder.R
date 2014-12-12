## peakfinder.R ##
#	Author: APM
# Date: 15Aug2013
#	Purpose:		utility functions for creating genetic maps from (interval) event data, and
#							for finding hot/cold spots on those maps.
# Thanks to:	Gary Churchill's scripts for G2F1 manuscript (Liu & Morgan (2013) PLoS Genetics), from which
#							most of this code was ported.

# external packages
library(GenomicRanges)
library(plyr)

## read.ranges() ##
# read a list of ranges specified samtools-style, like {chr}:{start}-{end}
read.regions <- function(file, ...) {
	df <- read.table(file, header = FALSE, sep = ":")
	df <- cbind( df, apply(simplify2array( strsplit(as.character(df[,2]), "-") ), 1, as.numeric ) )
	df[,2] <- NULL
	colnames(df) <- c("seqnames","start","end")
	return(df)
}

## integrate.events() ##
# get cumulative distribution of events along chromosomes, treating events as INTERVALS and not points
# allows for overlapping events; assumes intervals are HALF-OPEN by default (but this can be overridden)
# <events> = GRanges of event intervals, or a dataframe coercible to same
# <half.open> = assume intervals are half-open; if FALSE, they will be trimmed by 1 from the left
# integrate.events <- function(events, half.open = TRUE, weight.col = NULL, ...) {
# 	
# 	weighted <- !is.null(weight.col)
# 	
# 	if (is.data.frame(events)) {
# 		# handle dataframe input gracefully, provided it is GRanges()-ish
# 		gr <- with(events, GRanges(seqnames = Rle(seqnames),
# 															 ranges = IRanges(start = start, end = end)) )
# 		if (weighted) {
# 			values(gr)[ ,weight.col ] <- events[ ,weight.col ]
# 		}
# 		events <- gr
# 	}
# 	
# 	if (!half.open) {
# 		start(events) <- start(events) + 1
# 	}
# 	
# 	cat("Starting integration...\n")
# 	
# 	ldply(seqlevels(events), function(chr) {
# 		
# 		gr <- events[ seqnames(events) == chr ]
# 		if ( !length(gr) ) return(NULL)
# 		
# 		if ( any(width(gr) == 0) ) {
# 			stop("Interval widths not strictly positive: bad news.")
# 		}
# 		
# 		pos <- sort(unique( c(start(gr), end(gr)) ))
# 		# unique.gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos[-length(pos)], end = pos[-1]))
# 		# print(pos)
# 		tally <- rep(0, length(pos))
# 		for (i in 2:length(pos)) {
# 			weights <- 1
# 			right <- pos[i]
# 			left <- pos[i-1]
# 			overlapping.events <- gr[ start(gr) <= left & end(gr) >= right ]
# 			denom <- (width(overlapping.events)-1)
# 			if (weighted) {
# 				weights <- values(overlapping.events)[ ,weight.col ]
# 			}
# 			delta <- sum( (right - left)*weights/(denom) )
# 			tally[i] <- tally[i-1] + max(0, delta)
# 			# print(paste(left, right))
# 			# print(paste("Width of overlapping intervals:", paste(denom, collapse = ",")))
# 			# if (length(overlapping.events) > 1) print(overlapping.events)
# 			# print(paste("Additional density:", delta))
# 			# print(paste("Running total:", tally[i]))
# 			# print("----------")
# 		}
# 		
# 		return( data.frame(seqnames = chr, pos = pos, cumsum = tally) )
# 		
# 	}, .progress = "text")
# 
# }

integrate.events <- function(events, ...) {
	
	strand(events) <- "*"
	bins <- disjoin(events)
	print(bins)
	ol <- countOverlaps(bins, events)
	weight <- width(bins)*ol
	#print(head(ol, 100))
	#print(length(events))
	#print(length(bins))
	
	values(bins)$tally <- weight
	btwn <- gaps(bins)
	values(btwn)$tally <- 0
	rez <- sort(c(btwn, bins))
	values(rez)$cumsum <- cumsum(values(rez)$tally)
	return(rez[ strand(rez) == "*" & as.vector(seqnames(rez)) %in% as.vector(seqnames(events)) ])
	
}

## pad.map() ##
# utility for padding a cumulative map (eg. from integrate.events()) so that it covers appropriate region of chromosome
# for recombination maps, this region should be the part of the chromosome between the first and last marker per chrom on the array
# <map> = dataframe with cumulative map (eg. from integrage.events()) which is assumed to be sorted in a sensible order
# <starts> = named vector of leftmost positions on each chromosome
#	<ends> = named vector of rightmost positions on each chromosome
# intended for use with ddply:
#		ddply(map, .(seqnames), pad.map, seqlengths)
pad.map <- function(map, starts, ends, ...) {
	
	# sanity check
	stopifnot( length(starts) == length(ends) & all(ends >= starts) )
	
	first.row <- map[1,]
	first.row$pos <- starts[ as.character(map$seqnames[1]) ]
	last.row <- map[ nrow(map), ]
	last.row$pos <- ends[ as.character(map$seqnames[1]) ]
	
	return( rbind(first.row, map, last.row) )
	
}

## .interpolateFromBp() ##
# [internal use only]
#	utility for linear interpolation of numeric values within integer intervals along genome
# originally by Gary Churchill 10Jan2012 for interpolating map distances (cM), but is general
# <Map> = dataframe with 2 columns, first is position and second is value to be interpolated
.chrInterpolateFromBp <- function(Map, bpPositions) {
	foundIndices <- findInterval(bpPositions, Map[, 1])
	
	interpPositions <- c()
	mapRowCount <- nrow(Map)
	for(posIndex in 1 : length(bpPositions))
	{
		pos <- bpPositions[posIndex]
		foundIndex <- foundIndices[posIndex]
		if(foundIndex == 0) {
			interpPositions[posIndex] <- NA
		} else {
			foundBpPos <- Map[foundIndex, 1]
			foundCmPos <- Map[foundIndex, 2]
			if(pos == foundBpPos) {
				interpPositions[posIndex] <- foundCmPos
			} else {
				nextBpPos <- Map[foundIndex + 1, 1]
				nextCmPos <- Map[foundIndex + 1, 2]
				bpGap <- nextBpPos - foundBpPos
				cmGap <- nextCmPos - foundCmPos
				
				gapWeight <- (pos - foundBpPos) / bpGap
				interpPositions[posIndex] <- foundCmPos + (cmGap * gapWeight)
			}
		}
	}
	interpPositions
}

## make.density() ##
#	calculate event density from cumulative distribution using sliding window
# <map> = dataframe from integrate.events(); expected to have columns "seqnames" and "pos" and to be sorted
# intended for use with ddply:
#		ddply(map, .(seqnames), make.density)
make.density <- function(map, count.col = "cumsum", step = 100e3, window = 1000e3, ...) {
	
	if (!is.data.frame(map)) {
		map <- as.data.frame(map)
		map$pos <- map$start
	}

	rez <- NULL
	if (nrow(map) > 1) {
		
		start <- seq(min(map$pos), max(map$pos) - window, step)
		end <- start + window
		# print(cbind(start, end))
	
		events.start <- .chrInterpolateFromBp(map[ ,c("pos",count.col) ], start)
		events.end <- .chrInterpolateFromBp(map[ ,c("pos",count.col) ], end)
		
		rez <- data.frame(seqnames = map$seqnames[1], pos = (start + end)/2, density = (events.end - events.start)/window )
		
	}
	
	return(rez)
	
}

## .forward.sw() ##
# [internal use only]
# recursion function to do a 1-D Smith-Waterman (1d-SW) forward-pass on precomputed enrichment scores
# and returns a dataframe with new element "E"
# originally by Gary Churchill 10Jan2012
# <map> = dataframe sorted in a meaningful order
# <score.col> = column in <map> which contains the score on which to do 1d-SW
.forward.sw <- function(map, score.col = "e.score", ...) {
	
	scores <- as.numeric(map[ ,score.col ])
	E <- rep(0, nrow(map))
	E[1] <- max(0, scores)
	
	for(i in 2:length(scores)){
		E[i] <- max(0, E[i-1] + scores[i])
	}
	
	map$E <- E
	return(map)
	
}

## make.excursion() ##
# compute an "excursion score", the result of a 1d-SW pass on enrichment scores
# <map> = dataframe with columns (at least) "seqnames", "pos", <density.col>
# <theta> = enrichment factor for hotspot (eg. 100) or coldspot (eg. 1/100)
# <density.col> = column name which holds the value on which enrichment score is computed
#
# return value is a new dataframe with added columns:
# <e.score> = enrichment score for each window
# <E> = result of Smith-Waterman forward pass = excursion score
# designed for use with ddply by chromosome:
#		ddply(map, .(seqnames), 1000, 1/100)
make.excursion <- function(map, density.col = "density", theta = 100) {
	
	# compute the e.scores
	dens <- map[ ,density.col ]
	lambda <- mean(dens)
	gamma <- lambda * theta
	map <- transform( map, e.score = lambda - gamma + dens * log(gamma/lambda) )
	
	# compute the forward pass
	map <- .forward.sw(map, "e.score")
	return(map)

}

# test it
make.excursion(df2)

## get.peaks() ##
# extract best segment(s) from an excursion plot
# taking a list of scores and the output of the 1-D forward pass, ie. the output from make.excursion(),
# it extracts the best segment and resets the scores over that segment to zero.
# Following extraction of the best segment, the forward pass can be re-executed 
# and the next segment extracted...
#
# <map> = dataframe sorted in sensible order, with columns "seqnames", "pos", <score.col> (enrichment score) and "E" (excursion score)
# <minscore> = score cutoff below which function will quit trying to find peaks (not sure how to interpret???)
# <score.col> = column name holding enrichment scores (for use with .forward.sw())
#
# return value is a dataframe suitable for conversion to a GRanges, with columns "seqnames","start","end","E"
get.peaks <- function(map, minscore = 10, score.col = "e.score") {
	
	peaks <- NULL
	E <- map$E
	while(max(E) > minscore) {
		
		# find index of the maximum excursion
		E <- map$E
		idx.end <- max(which(E == max(E)))

		# find the last zero preceeding the max index
		zero.list <- which(E[1:idx.end] == 0)
		idx.start <- ifelse( length(zero.list) == 0, 1, max(zero.list) )
		
		# report region
		tmp <- cbind( map[ idx.start,c("seqnames","pos") ], map[ idx.end,c("pos","E") ] )
		names(tmp) <- c("seqnames","start","end","E")
		peaks <- rbind(peaks, tmp)

		#reset the scores to zero
		map[ idx.start:idx.end,score.col ] <- 0

		print("found a peak; going again...")
		# recompute the forward pass
		map <- .forward.sw(map, score.col)
	}
	
	# sort peaks by genome coordinate
	if (!is.null(peaks)) {
		if(nrow(peaks) > 1) {
			peaks <- peaks[ order(peaks$start), ]
		}
	}
	
	return(peaks)

}