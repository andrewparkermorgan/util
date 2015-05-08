## miscellaneous functions for plotting in genomic coordinates with ggplot2

library(GenomicRanges)
library(ggplot2)
library(scales)
library(grid)

## scales & transformations
trans_genome <- function(scale = 1e6, ...) {
	trans_new("genome coordinates", transform = function(x) x/scale, 
						inverse = function(X) X * scale)
}

scale_x_genome <- function(scale = 1e6, unit = NULL, ...) {
	
	if (is.null(unit)) {
		if (scale == 1e3)
			unit <- "kbp"
		else if (scale == 1e6)
			unit <- "Mbp"
		else if (scale == 1e9)
			unit <- "Gbp"
	}
	
	gtrans <- trans_genome(scale, unit)
	return(list( scale_x_continuous(trans = gtrans, labels = function(x) gtrans$transform(x), ...),
							 xlab(paste0("\nposition (", unit, ")")) ))

}

neglog_trans <- function (base = exp(1)) 
{
	trans <- function(x) -1*log(x, base)
	inv <- function(x) base^(-1*x)
	trans_new(paste0("-log-", format(base)), trans, inv, log_breaks(base = base), 
			  domain = c(1e-100, Inf))
}
	

scale_y_logp <- function(...) {
	
	scale_y_continuous(..., trans = neglog_trans(10))
	
}

## themes
theme_nothing <- function(base_size = 12, base_family = "Helvetica")
{
	theme_bw(base_size = base_size, base_family = base_family) %+replace%
		theme(
			rect             = element_blank(),
			line             = element_blank(),
			text             = element_blank(),
			axis.ticks.margin = unit(0, "lines")
		)
}

theme_gbrowse <- function(base_size = 12, base_family = "Helvetica") {
	theme_bw(base_size = base_size, base_family = base_family) %+replace%
		theme(
			panel.grid = element_blank(),
			panel.border = element_blank(),
			axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			axis.ticks.margin = unit(0.25, "lines")
		)
}

theme_axesonly <- theme_gbrowse <- function(base_size = 12, base_family = "Helvetica") {
	theme_bw(base_size = base_size, base_family = base_family) %+replace%
		theme(
			panel.grid = element_blank(),
			panel.border = element_blank(),
			axis.line.y = element_blank(),
			#axis.ticks.y = element_blank(),
			#axis.text.y = element_blank(),
			#axis.title.y = element_blank(),
			axis.ticks.margin = unit(0.25, "lines")
		)
}

## S3 method overload to autoplot a GRanges with sensible axes
## geom is left up to the user
ggplot.GRanges <- function(gr, vspace = 0.2, geom = "rect", ....) {
	
	df <- as.data.frame(gr)
	df$seqnames <- factor(df$seqnames, levels = rev(levels(df$seqnames)))
	df$z <- as.numeric(df$seqnames)
	df$vspace <- vspace
	
	p <- ggplot(df, aes(xmin = start, xmax = end, ymin = z-(1-vspace/2), ymax = z-(vspace/2))) +
		scale_y_continuous(breaks = seq_along(levels(df$seqnames))-0.5, labels = levels(df$seqnames)) +
		scale_x_genome() +
		theme_bw()
	if (geom == "rect") {
		p <- p + geom_rect()
	}
	
	return(p)
	
}

## geoms
geom_tx <- function(exons, at = 0, height = 1, arrows = grid::arrow(length = unit(4, "pt"), type = "closed"), do.introns = TRUE, ...) {
	
	if (!inherits(exons, "GRanges"))
		exons <- makeGRangesFromDataFrame(exons, keep.extra.columns = TRUE)
	
	## get exons from introns
	.introns <- gaps(exons)[-1]
	introns <- droplevels(as.data.frame(.introns))
	
	if (any(strand(exons) == "-")) {
		x <- introns$start
		introns$start <- introns$end
		introns$end <- x
	}
	
	values(exons)$at <- at
	values(exons)$height <- height
	introns$at <- at
	introns$height <- height
	
	#print(exons)
	#print(introns)
	
	list( geom_rect(data = droplevels(as.data.frame(exons)), aes(xmin = start, xmax = end, ymin = at, ymax = at+height), ...),
				geom_segment(data = introns, aes(x = start, xend = end, y = at+height/2, yend = at+height/2), arrow = arrows, ...) )
	
}

## Haploview-style local LD heatmap using ggplot2
## <ld> = dataframe of pairwise LD as computed by plink --r2
gghaploview <- function(ld, ...) {
	
	x <- sort(unique(c(ld$BP_A, ld$BP_B)))
	ld$x <- match(ld$BP_A, x)
	ld$y <- match(ld$BP_B, x)
	#R <- sparseMatrix(i = ld$x, j = ld$y, x = ld$R2, giveCsparse = FALSE, index1 = TRUE)
	R <- matrix(0.01, nrow = length(x), ncol = length(x))
	R[ as.matrix(ld[ ,c("x","y") ]) ] <- ld$R2
	#ii <- R@i
	#jj <- R@j
	nz <- which(R > 0, arr.ind = TRUE)
	ii <- nz[,1]
	jj <- nz[,2]
	
	message("Matrixifying LD...")
	polys <- data.frame(x = numeric(4*length(ii)), y = numeric(4*length(ii)),
						id = numeric(4*length(ii)), r2 = numeric(4*length(ii)))
	pb <- txtProgressBar(min = 0, max = length(ii), style = 3)
	for (z in 1:length(ii)) {
		
		if (ii[z] >= jj[z])
			next
		
		polys$id[ 4*(z-1) + 1:4 ] <- z
		polys$r2[ 4*(z-1) + 1:4 ] <- R[ ii[z],jj[z] ]
		#polys$r2[ 4*(z-1) + 1:4 ] <- ld$R2[z]
		
		polys$y[ 4*(z-1)+1 ] <- ii[z]
		polys$x[ 4*(z-1)+1 ] <- jj[z]-1
		
		polys$y[ 4*(z-1)+2 ] <- ii[z]
		polys$x[ 4*(z-1)+2 ] <- jj[z]
		
		polys$y[ 4*(z-1)+3 ] <- ii[z]+1
		polys$x[ 4*(z-1)+3 ] <- jj[z]
		
		polys$y[ 4*(z-1)+4 ] <- ii[z]+1
		polys$x[ 4*(z-1)+4 ] <- jj[z]-1
		
		setTxtProgressBar(pb, z)
		
	}
	
	theta <- -pi/4
	rot <- matrix(c(cos(theta), -sin(theta),
					sin(theta), cos(theta)),
					ncol = 2, byrow = TRUE)
	polys[ ,1:2 ] <- t( rot %*% t(polys[,1:2 ]) )
	bounds <- range(x)
	ggplot(polys) +
		geom_polygon(aes(x = x, y = y-2, group = id, fill = r2), colour = NA) +
		annotate("point", pch = 19, x = min(polys$x) + max(polys$x)*(x-min(bounds))/diff(bounds),
				 y = 0.02*max(polys$x)) +
		annotate("segment", x = (1:length(x)*sqrt(2))-0.0*sqrt(2),
				 xend = min(polys$x) + max(polys$x)*(x-min(bounds))/diff(bounds),
				 y = -1, yend = 0.02*max(polys$x)) +
		geom_hline(yintercept = 0.02*max(polys$x)) +
		annotate("segment", x = max(polys$x)*(pretty(x)-min(bounds))/diff(bounds),
				 xend = max(polys$x)*(pretty(x)-min(bounds))/diff(bounds),
				 y = 0.02*max(polys$x), yend = 0.05*max(polys$x)) +
		annotate("text", label = pretty(x)/1e6, x = max(polys$x)*(pretty(x)-min(bounds))/diff(bounds), y = 0.1*max(polys$x)) +
		scale_fill_gradient(expression(widehat(r^2)), low = "lightyellow", high = ("red")) +
		coord_equal() +
		theme_gbrowse() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
	
}

gghaploview2 <- function(ld, min.r2 = 0.1, colours = c("lightyellow","yed"), gridlines = NA, ...) {
	
	require(data.table)
	if (!is.data.table(ld))
		ld <- data.table(ld)
	
	mk.a <- unique( ld[ ,.(SNP_A, BP_A) ] )
	setnames(mk.a, c("marker","pos"))
	mk.b <- unique( ld[ ,.(SNP_B, BP_B) ] )
	setnames(mk.b, c("marker","pos"))
	mk <- unique(rbind(mk.a, mk.b))
	mk$o <- order(mk[,pos])-1
	setkeyv(mk, "marker")
	
	message("Creating LD matrix...")
	allpolys <- data.table()
	corners <- matrix( c(0,1,
						 0,0,
						 1,0,
						 1,1),
						 byrow = TRUE, ncol = 2 )
	for (i in 1:nrow(corners)) {
		ii <- as.character(ld$SNP_A)
		jj <- as.character(ld$SNP_B)
		x <- as.numeric(mk[ ii,o ]) + corners[i,1]
		y <- as.numeric(mk[ jj,o ]) + corners[i,2]
		r2 <- ld$R2
		id <- 1:nrow(ld)
		polys <- data.table(x = x, y = y, r2 = r2, id = id, marker = jj, corner = i)
		allpolys <- rbind(allpolys, polys)
	}
	
	message("Rotating LD matrix...")
	theta <- -pi/4
	rot <- matrix(c(cos(theta), -sin(theta),
					sin(theta), cos(theta)),
					ncol = 2, byrow = TRUE)
	X <- as.matrix(allpolys[ ,.(x,y) ])
	Xt <- t( rot %*% t(X) )
	allpolys[ ,x := as.vector(Xt[,1]) ]
	allpolys[ ,y := as.vector(Xt[,2]) ]
	setkeyv(allpolys, c("id","corner"))
	#return(allpolys)
	
	pos <- sort(as.numeric(mk[allpolys[ ,min(x), by = "marker" ],pos]))
	diagpos <- allpolys[ ,min(x), by = "marker" ][,V1]
	bounds <- range(pos)
	maxpos <- allpolys[ ,max(x) ]
	minpos <- allpolys[ ,min(x) ]
	SCALE <- 0.1
	baseline <- minpos-SCALE*maxpos
	
	message("Drawing plot...")
	ggplot(as.data.frame(allpolys[ r2 > min.r2 ])) +
		# draw LD squares
		geom_polygon(aes(group = id, x = x, y = y, fill = r2), colour = gridlines) +
		# draw markers
		annotate("point", pch = 19, x = minpos + maxpos*(pos-min(bounds))/diff(bounds),
				 y = baseline) +
		# draw lines between markers and LD squares
		annotate("segment", xend = diagpos*2,
				 x = minpos + maxpos*(pos-min(bounds))/diff(bounds),
				 y = baseline, yend = baseline + 0.75*SCALE*maxpos) +
		# draw horizontal axis
		geom_hline(yintercept = baseline) +
		# draw axis ticks
		annotate("segment", x = maxpos*(pretty(pos)-min(bounds))/diff(bounds),
				 xend = maxpos*(pretty(pos)-min(bounds))/diff(bounds),
				 y = baseline, yend = baseline - 0.25*SCALE*maxpos) +
		# draw axis labels
		annotate("text", label = pretty(pos)/1e6, x = maxpos*(pretty(pos)-min(bounds))/diff(bounds),
				 y =  baseline - 0.5*SCALE*maxpos) +
		scale_fill_gradientn(expression(widehat(italic(r)^2)), colours = colours) +
		coord_equal() +
		theme_gbrowse() + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
	
}

ggmanhattan <- function(df, seqnames = NULL, space = 10, cols = c("grey30","grey60"), ...) {
	
	if (all(c("chr","pos") %in% colnames(df)))
		df <- df[ ,c("chr","pos", setdiff(colnames(df), c("chr","pos"))) ]
	
	if (!is.factor(df[,1]))
		if (!is.null(seqnames))
			df[,1] <- factor(df[,1], levels = seqnames)
		else
			df[,1] <- factor(df[,1])
	else
		df[,1] <- factor(df[,1])
	
	chrlen <- tapply(df[,2], df[,1], max, na.rm = TRUE)/1e6 + space
	adj <- c(0, cumsum(chrlen))
	ticks <- adj[ -length(adj) ] + diff(adj)/2
	#names(adj) <- c(names(chrlen),"z")
	#print( cbind(chrlen, adj[-1], ticks) )
	#print(length(chrlen))
	#print(length(ticks))
	df$.adj <- adj[ as.numeric(df[,1]) ]
	df$.x <- df$.adj + df[,2]/1e6
	df$.chr <- df[,1]
	colmap <- setNames( rep_len(c(0,1), nlevels(df$.chr)), levels(df$.chr)  )
	df$.colour <- factor(colmap[ df$.chr ])
	#print(tapply(df$.x, df$.chr, max))
	
	rez <- ggplot(df, aes(x = .x, colour = .colour))
	rez <- rez +
		scale_x_continuous(breaks = ticks, minor_breaks = adj,
						   labels = gsub("^chr", "", names(chrlen))) +
		scale_colour_manual(values = cols) +
		guides(colour = FALSE) +
		theme_bw() + theme(axis.title.x = element_blank(),
						   panel.grid.minor = element_line(colour = "grey90"),
						   panel.grid.minor = element_blank())
	return(rez)
	
}

plot.haplotypes <- function(haps, flatten = TRUE, fill = "strain", omit.chrs = c("chrY","chrM"), shorten.chrs = FALSE, ...) {
	
	## flatten list-of-lists, if necessary
	haps.df <- data.frame()
	if (flatten | !is.data.frame(haps))
		haps.df <- ldply(haps, as.data.frame)
	else
		haps.df <- haps
	
	## keep only desired chromosomes (1-19,X by default)
	haps.df <- subset(haps.df, !(seqnames %in% omit.chrs))
	chroms <- setdiff( levels(haps.df$seqnames), omit.chrs )
	if (shorten.chrs) {
		chroms <- gsub("^chr","",chroms)
		haps.df$seqnames <- gsub("^chr","",as.character(haps.df$seqnames))
	}
	haps.df$seqnames <- factor( haps.df$seqnames, levels = chroms )
	
	haps.df$seqnames <- with(haps.df, factor(as.character(seqnames), levels = rev(levels(seqnames))) )
	haps.df$updown <- ifelse(haps.df$origin == levels(haps.df$origin)[1], 1, -1)
	haps.df$ypos <- as.numeric(haps.df$seqnames)
	haps.df <- transform(haps.df, ymin = 0.4*updown+ypos, ymax = 0.05*updown+ypos)
	haps.df$fill <- haps.df[ ,fill ]
	
	p <- ggplot(haps.df) +
		geom_rect(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = fill)) +
		scale_y_continuous(breaks = 1:length(chroms), labels = rev(chroms)) +
		scale_x_genome() +
		theme(axis.title.y = element_blank(), plot.title = element_text(face = "bold"))
	
	if ("id" %in% colnames(haps.df))
		p <- p + ggtitle(paste0(haps.df$id[1], "\n"))
	
	return(p)
	
}

warped.breaks <- function(x, y1, y2, scale = c("native","transformed"), ...) {
	
	stopifnot(length(y1) == length(y2))
	i <- findInterval(x, y1)
	b <- pretty(y2[i])
	b <- b[ b > min(y2) & b < max(y2) ]
	
	scale <- match.arg(scale)
	if (scale == "native")
		return(y1[ findInterval(b, y2) ])
	else
		return(b)
}

## plot genotypes (from genotypes object) as evenly-spaced circles, with actual positions indicated on x-axis below
dotplot.genotypes <- function(gty, size = 2, meta = NULL, ...) {
	
	require(RColorBrewer)
	require(plyr)
	require(reshape2)
	
	## check that input is right sort
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes.'")
	
	if (is.null(attr(gty, "map")))
		stop("Please attach a marker map to the input object.")
	
	## set up horizontal positions of markers
	map <- attr(gty, "map")
	map$i <- as.numeric(factor(map$pos))
	map$ipos <- map$pos - min(map$pos)
	map$relpos <- map$ipos/diff(range(map$pos)) * max(map$i)
	map <- map[ with(map, order(i)), ] # important
	#print(head(map))
	
	## convert genotypes to long-form for ggplot
	df <- melt(gty)
	colnames(df) <- c("marker","id","ibs")
	if (is.numeric(gty))
		df$ibs.code <- factor(df$ibs)
	else
		df$ibs.code <- factor(df$ibs, levels = c("A","C","G","T","H"))
	df <- merge(df, map)
	
	## add family info, if present
	if (!is.null(attr(gty, "ped")))
		df <- merge(df, attr(gty, "ped"), by.x = "id", by.y = "iid", all.x = TRUE)
	
	## add metadata, if provided
	if (!is.null(meta))
		df <- merge(df, meta, all.x = TRUE)
	
	## check that dimensions still match
	stopifnot(nrow(df) == prod(dim(gty)))
	print(head(df))
	
	p <- ggplot(subset(df, !is.na(ibs))) +
		geom_tile(aes(x = i, y = id, fill = ibs.code), pch = 21, size = size) +
		geom_segment(data = map, aes(x = i, xend = relpos, y = 0.8, yend = 0), colour = "grey70") +
		geom_hline(yintercept = 0) +
		scale_x_continuous(breaks =  warped.breaks(map$relpos, map$relpos, map$pos),
						   labels = round(warped.breaks(map$relpos, map$relpos, map$pos, scale = "transformed")/1e6, 2)) +
		theme_minimal() + theme(axis.title.y = element_blank()) +
		xlab("\nposition (Mbp)")
	
	if (is.numeric(gty))
		p <- p + scale_fill_manual("IBS", values = c("white","grey","black"))
	else
		p <- p + scale_fill_manual("genotype", values = brewer.pal("Set1", 9)[ c(1:4,9) ])
	
	attr(p, "map") <- map
	
	return(p)
	
}