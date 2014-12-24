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

## themes
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

## geoms
geom_tx <- function(exons, at = 0, height = 1, arrows = grid::arrow(length = unit(4, "pt"), type = "closed"), ...) {
	
	if (!inherits(exons, "GRanges"))
		exons <- makeGRangesFromDataFrame(exons, keep.extra.columns = TRUE)
	
	## get exons from introns
	.introns <- gaps(exons)[-1]
	introns <- as.data.frame(.introns)
	
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
	
	list( geom_rect(data = as.data.frame(exons), aes(xmin = start, xmax = end, ymin = at, ymax = at+height), ...),
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
	
	polys <- data.frame(x = numeric(4*length(ii)), y = numeric(4*length(ii)),
						id = numeric(4*length(ii)), r2 = numeric(4*length(ii)))
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