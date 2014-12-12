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
