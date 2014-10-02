plot.stacked.haps <- function(.haps, chroms = "chr2", which = "both",
															sort.col = 2, sort.on = "recombinant", manual.sort = NULL,
															strain.col = "strain", label = TRUE, xlim = NULL, ...) {
	
	require(ggplot2)
	require(plyr)
	require(GenomicRanges)
	
	## copy original input
	haps <- .haps
	
	if (!inherits(haps, "haplotypes")) {
		warning( paste0("Plotting method assumes data structure provided for class 'haplotypes' defined by Rmegamuga package",
										"(although the package itself isn't require()d by this method)") )
	}
	
	## set target range (for clipping segments which fall outside plotting limits)
	x.target <- GRanges()
	if (!is.null(xlim))
		x.target <- GRanges(seqnames = chroms, ranges = IRanges(start = xlim[1], end = xlim[2]))
	
	if (length(which) == length(haps)) {
		## given a vector <which> parallel to <haps>, use its i-th element to subset the i-th segment set
		haps <- lapply(1:length(which), function(i) {
			rez <- haps[[i]][[ which[i] ]]
			if (!is.null(rez))
				rez[ rez %over% x.target ]
			else
				return(x.target)
		})
		names(haps) <- names(.haps)
	}
	else if (which[1] == "both") {
		## mode "both": now figure out which haplotype (mat/pat) to use for sorting
		haps <- lapply(haps, function(h) {
			# first limit segments to the desired chromosome and interval
			h <- lapply(h, function(this.h) {
				if (is.null(xlim)) {
					xlim <- c(1, seqlengths(this.h)[ chroms ])
					x.target <- GRanges(seqnames = chroms, ranges = IRanges(start = xlim[1], end = xlim[2]))
				}
				this.h <- this.h[ this.h %over% x.target ]
				return(this.h)
			})
			# now decide which to keep (the first in the list, by default)
			i <- 1
			if (sort.on == "recombinant") {
				# mode "recombinant": sort using the one with the MOST recombinatinos
				i <- which.max(sapply(h, length))
			}
			else if (sort.on == "first") {
				# mode "first": sort using the one with the FIRST (leftmost) recombination
				first <- sapply(h, function(this.h) {
					starts <- start(h)
					min(starts[ starts > xlim[1] & starts < xlim[2] ])
				})
				i <- max(which.min(first), 1)
			}
			rez <- unlist(GRangesList(h[ i[1] ]))
			# values(rez)$origin <- names(h)[ i[1] ]
			return( rez )
		})
	}
	else {
		## subset using first element of <which>
		haps <- lapply(haps, "[[", which[1])
	}
	haps <- haps[ !sapply(haps, function(h) is.null(h) | !length(h)) ]
	# print(haps)
	# print(sapply(haps, length))
	
	## sort chromosomes
	sort.by <- matrix(1:length(haps), nrow = length(haps), ncol = 1)
	rownames(sort.by) <- names(haps)

	if (!is.null(manual.sort)) {
		sort.by <- matrix(manual.sort, nrow = length(manual.sort), ncol = 1)
		rownames(sort.by) <- names(manual.sort)
		sort.col <- 1
	}
	else if (sort.col > 0) {
		
		# sort by recombination breakpoint(s)
		# do it recursively, to handle multiply-recombinant chromosomes: by 0th ... nth event
		n.recombs <- sapply(haps, length) - 1
		print(n.recombs)
		last.recomb <- 0
		pos <- matrix(0, nrow = length(haps), ncol = max(n.recombs)+1)
		while( any(n.recombs >= last.recomb) ) {
			i <- which(n.recombs >= last.recomb)
			pos[ i,(last.recomb+1) ] <- sapply(haps[i], function(h) start(h)[ (last.recomb+1) ])
			last.recomb <- last.recomb + 1
			# print(pos)
		}
		sort.by <- cbind(sort.by, pos)
		o <- order(sort.by[ ,sort.col ])
		sort.by <- sort.by[o,]
	}
	# print(sort.by)
	sort.by[,1] <- 1:nrow(sort.by)
	final.sort.order <- sort.by[,1]
	# print(final.sort.order)
	
	## roll up haplotypes into dataframe for plotting, adding sort order
	haps.df <- data.frame()
	if (length(which) == 1) {
		attr(.haps, "split_labels") <- data.frame(id = names(.haps)) # coax plyr into adding ids automatically...
		haps.df <- ldply(.haps, function(h) ldply(h, as.data.frame))
		colnames(haps.df)[1:2] <- c("id","origin")
	}
	else {
		haps.df <- ldply(haps, as.data.frame)
		haps.df$origin <- which[ as.character(haps.df$.id) ]
		colnames(haps.df)[ c(1, ncol(haps.df)) ] <- c("id","origin")
	}
	
	print(head(haps.df))
	if (!is.null(xlim)) {
		## clip segments which have an endpoint in the target interval
		haps.df$start <- pmax(haps.df$start, xlim[1])
		haps.df$end <- pmin(haps.df$end, xlim[2])
	}
	haps.df$ypos <- final.sort.order[ as.character(haps.df$id) ]
	if (which == "both") {
		haps.df$updown <- ifelse(haps.df$origin == "maternal", 1, -1)
	}
	else {
		if (!is.null(xlim))
			# haps.df <- subset(haps.df, (start > xlim[1] & start < xlim[2]) |
			#										(end > xlim[1] & end < xlim[2]) )
		# haps.df <- subset(haps.df, origin == which)
		haps.df$updown <- 1
		haps.df$ypos <- haps.df$ypos - 0.5
	}
	spacing <- 0.1
	shrink <- ifelse(which[1] == "both", 0.5, 1)
	haps.df$ymin <- with(haps.df, ypos + pmin(0, updown*(shrink-spacing/2)))
	haps.df$ymax <- with(haps.df, ypos + pmax(0, updown*(shrink-spacing/2)))
	print(head(haps.df))
	
	## construct stacked-haplotype plot
	rez <- ggplot(haps.df) +
		geom_rect(aes(xmin = start/1e6, xmax = end/1e6, ymin = ymin, ymax = ymax, fill = strain)) +
		xlab("\nposition (Mbp)") +
		scale_y_continuous(breaks = final.sort.order, limits = c(0,max(final.sort.order)+1)) +
		theme_bw() +
		theme(axis.title.y = element_blank())
	if (!label) {
		rez <- rez + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
	}
	
	## tack on the sorting order, for use outside this function
	attr(rez, "sort.by") <- final.sort.order
	
	return(rez)
	
}

plot.recombinations <- function(recombs, chroms = "chr2", which = "paternal", target = NULL, pheno = NULL, pheno.col = 2, ...) {
	
	require(plyr)
	require(ggplot2)
	require(GenomicRanges)
	require(gridExtra)
	
	## get just the requested haplotype and chromosome(s)
	recombs <- lapply(recombs, "[[", which)
	recombs <- lapply(recombs, function(r) {
		if (!is.null(target)) {
			r[ r %over% target ]
		} else {
			r[ seqnames(r) %in% chroms ]
		}
	})
	
	recombs.df <- ldply(recombs, as.data.frame)
	recombs.df$id <- with(recombs.df, reorder(.id, start, min))
	
	## merge in phenotype data, if provided
	# NB: it must have a column 'id' for indexing
	if (!is.null(pheno)) {
		pheno$pheno.col <- pheno[,pheno.col]
		recombs.df <- merge(recombs.df, pheno, by = "id")
	}
	
	# recombs.df$stack <- base::order(recombs.df$start)
	# print(recombs.df)
	# recombs.df$group <- with(recombs.df, from:to)
	print(head(recombs.df))
	
	p1 <- ggplot(recombs.df) +
		geom_point(aes(x = (end+start)/2e6, y = id, colour = from:to)) +
		xlab("\nposition (Mbp)") +
		guides(fill = guide_legend("strains at junction")) +
		theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
	
	if (!is.null(pheno)) {
		
		p2 <- ggplot(recombs.df) +
			geom_point(aes(x = pheno.col, y = id)) +
			xlab(paste0("\n", pheno.col)) +
			theme_bw() + theme(axis.title.y = element_blank())
		return(list(recombs = p1, pheno = p2, data = recombs.df))
		
	}
	else {
		return(list(recombs = p1, data = recombs.df))
	}
	
	# ggplot(recombs.df) +
	# 	stat_ecdf(aes(x = (start+end)/(2*1e6), colour = group)) +
	# 	ylab("cumulative proportion of recombination events\n") + xlab("\nposition (Mbp)") +
	# 	theme_bw()
}

overlap.with.meta <- function(query, subject, metavars = colnames(values(query)), ...) {
	
	require(GenomicRanges)
	
	flag <- (query %over% subject)
	meta.flag <- matrix(ncol = 0, nrow = length(query))
	
	if (length(subject) > 1 & !(length(subject) == length(query))) {
		warning("Result will probably not be well-defined due to shape of query vs subject ranges.")
	}
	
	for (i in colnames(values(query))) {
		if (i %in% colnames(values(subject))) {
			# hits.this.subject <- sapply(1:length(subject), function(j) flag & (as.character(values(query)[,i]) == as.character(values(subject)[j,i])) )
			# print(hits.this.subject)
			hits.this.subject <- flag & (values(query)[,i] == values(subject)[,i])
			meta.flag <- cbind(meta.flag, hits.this.subject)
		}
	}
	
	return(meta.flag)
	
}

roll.haps <- function(haps, ...) {
	
	require(GenomicRanges)
	require(plyr)
	
	dlply(haps, .(id), function(d) {
		rez <- dlply(d, .(origin), function(p) {
			with(p, GRanges(seqnames = Rle(seqnames), seqlengths = Rmegamuga:::seqlengths,
											ranges = IRanges(start = start, end = end),
											strain = strain, id = id))
		})
		class(rez) <- c("haplotypes", class(rez))
		return(rez)
	})
	
}

roll.recombinations <- function(haps, seqlengths = Rmegamuga:::seqlengths, ...) {
	
	require(GenomicRanges)
	require(plyr)
	
	forbidden <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
								 "isCircular", "start", "end", "width", "element")
	
	dlply(haps, .(id), function(d) {
		rez <- dlply(d, .(origin), function(p) {
			gr <- with(p, GRanges(seqnames = Rle(seqnames), seqlengths = seqlengths,
											ranges = IRanges(start = start, end = end),
											from = from, to = to, id = id))
			values(gr) <- cbind(values(gr), p[ ,setdiff(colnames(p), c(forbidden,colnames(values(gr)))) ])
			return(gr)
		})
		# class(rez) <- c("haplotypes", class(rez))
		return(rez)
	})
	
}

invert.recombinations <- function(recombs, seal = FALSE, ...) {
	
	require(GenomicRanges)
	
	forbidden <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
								 "isCircular", "start", "end", "width", "element")
	
	lapply(recombs, function(r.all) {
		segments <- GRanges()
		for (c in unique(runValue(seqnames(r.all)))) {
			r <- r.all[ seqnames(r.all) == c ]
			if (seal) {
				# 'seal' haplotype segments together -- eg. get rid of uncertainty around recomb events
				g <- GRanges(seqnames = c, seqinfo = seqinfo(r),
										 ranges = IRanges(start = c(1, pmin( mid(ranges(r))+1, seqlengths(r)[c] )),
										 								  end = c(mid(ranges(r)), seqlengths(r)[c])),
										 strain = c( as.character(values(r)$from), rev(as.character(values(r)$to))[1] )
										 )
				segments <- c(segments, g)
			}
			else {
				if (length(r)) {
					if (length(r) > 1) {
						i <- 1:(length(r)-1)
						# un-overlap possibly overlapping recomb events
						end(r)[i] <- pmin( end(r)[i], start(r)[i+1]-1 )
						start(r)[i+1] <- pmax( end(r)[i]+2, start(r)[i+1] )
					}
					r <- sort(r)
					g <- gaps(r)
					g <- g[ strand(g) == "*" & (as.vector(seqnames(g)) %in% as.vector(seqnames(r))) ]
					values(g)$strain <- c( as.character(values(r)$from), rev(as.character(values(r)$to))[1] )
					segments <- c(segments, g)
				}
			}
		}
		# print(segments)
		return(segments)
	})
	
}

other.haplotype <- function(x, coding = c("maternal", "paternal"), ...) {
	
	## flip haplotypes: if maternal, then return paternal, and vice versa
	
	which.hap <- outer(x, coding, "==")
	# print(which.hap)
	other.hap <- apply(which.hap, 2, "!")
	# print(other.hap)
	rez <-coding[ ((other.hap[,1]+1) %% 2) + 1 ]
	if (!is.null(names(x)))
		names(rez) <- names(x)
	
	return(rez)
	
}

## pseudo-phase DO-like haplotypes
pseudophase <- function(hap, hap1.col = "hap1", hap2.col = "hap2", seqlengths, ...) {
	
	require(GenomicRanges)
	
	swap <- function(row, j = 1, k = j + 1, ...) {
		row[j] <- tmp
		row[k] <- row[j]
		row[j] <- tmp
		return(row)
	}
	
	hap <- with(hap, hap[ order(start, end), ])
	if (nrow(hap) < 2) {
		for (i in 2:nrow(hap)) {
			this.hap1 <- hap[ i,hap1.col ]
			last.hap1 <- hap[ i-1,hap1.col ]
			this.hap2 <- hap[ i,hap2.col ]
			last.hap2 <- hap[ i-1,hap2.col ]
			if ((this.hap1 == last.hap2) | (this.hap2 == last.hap1)) {
				hap[i,] <- swap(hap[i,], hap1.col, hap2.col)
			}
		}
	}
	
	phased <- list(GRanges(), GRanges())
	names(phased) <- c(hap1.col, hap2.col)
	for (col in c(hap1.col, hap2.col)) {
		tmp <- Rle(hap[ ,col ])
		new.segments <- GRanges(seqnames = hap$seqnames[1], seqlengths = seqlengths,
														ranges = IRanges(start = hap[ start(tmp),"start" ], end = hap[ end(tmp),"end" ]),
														strain = runValue(tmp))
		phased[[col]] <- c(phased[[col]], new.segments)
	}
	
	return(phased)
	
}