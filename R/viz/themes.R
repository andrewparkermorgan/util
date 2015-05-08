## themes and miscellanious utility functions for ggplot2
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(RColorBrewer)

## Hadley's totally blank theme
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

## theme for horizontal dotplots, only y-axis ticks and text are present
theme_ylabs_only <- function(...) {
	
	thm <- theme_clean(...)
	thm + theme(axis.line = element_blank(), axis.title = element_blank(),
				axis.text.x = element_blank(), axis.ticks.x = element_blank() )
	
}

## a clean black-and-white theme with friendly way to specify legend placement
theme_clean <- function(legend.pos = c("outside","topright", "topleft", "bottomright", "bottomleft"), ...) {
	
	thm <- theme_classic(...)
	
	legend.pos <- match.arg(legend.pos)
	
	if (legend.pos == "outside") {
		thm <- thm + theme(legend.key.size = unit(0.8, "lines"))
	}
	else {
		lp <- lj <- c(0,0)
		if (legend.pos == "topleft") {
			lp <- lj <- c(0,1)
		}
		else if (legend.pos == "topright") {
			lp <- c(1,1)
			lj <- c(1,1)
		}
		else if (legend.pos == "bottomleft") {
			lp <- c(0,0)
			lj <- c(0,0)
		}
		else if (legend.pos == "bottomright") {
			lp <- lj <- c(1,0)
		}
		thm <- thm + theme(legend.position = lp, legend.justification = lj,
							legend.key.size = unit(0.8, "lines"),
							legend.background = element_blank() )
	}
	
	return( thm + theme(strip.background = element_blank(),
						strip.text = element_text(face = "bold")) )

}

## clean theme with x-axis only
theme_clean_xonly <- function(...) {
	
	theme_clean(...) %+replace% theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
									  axis.text.y = element_blank(), axis.title.y = element_blank())
	
}

## black-and-white theme with x-axis labels at 45 degree angle, like for a boxplot
theme_slanty_x <- function(base.theme = theme_bw, ...) {
	
	base.theme(...) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
						 axis.title.x = element_blank())
	
}

## clean theme for heatmaps
theme_heatmap <- function(...) {
	
	theme_classic(...) + theme(axis.text = element_text(size = 8),
							   axis.text.x = element_text(angle = 90, hjust = 1),
							   axis.line = element_blank(), axis.title = element_blank())
	
}

scale_fill_heatmap <- function(...) {
	
	require(RColorBrewer)
	
	heat.cols <- colorRampPalette(brewer.pal(4, "Spectral")[2:1])(4)
	scale_fill_gradientn(..., colours = heat.cols, na.value = "grey70")
	
}

## create a grob with a text label (like for labelling a plot panel),
##	suitable for placement with grid.arrange()
panel.label.grob <- function(label, pos = c("topleft","topright"),
							 y = 0.9, x = 0, hjust = 0, vjust = 1, ...) {
	
	pos <- match.arg(pos)
	if (pos == "topright") {
		x <- 0
		hjust <- 0
	}
	else if (pos == "topleft") {
		x <- 0.9
		hjust <- 1
	}
	
	grid.text(label, x = x, y = y, hjust = hjust, vjust = vjust, gp = gpar(...), draw = FALSE)
	
}

blank.grob <- function(...) {
	grid.rect(gp = gpar(col = NA))
}

## get an interpolator function for a RColorBrewer palette
brewer.interpolate <- function(palette, ...) {
	
	maxcol <- brewer.pal.info[ palette,"maxcolors" ]
	pal <- colorRampPalette(brewer.pal(maxcol, palette))
	return(pal)
	
}

## like scale_colour_brewer(), but allow arbitrary number colors via interpolation
scale_brewer_interpolated <- function(..., palette = "Spectral") {
	
	discrete_scale("colour", "brewer", palette = brewer.interpolate(palette), na.value = "grey50", ...)
	
}

## format number as pretty percentage
as.percent <- function(x) {
	percent_format()(x)
}

## arrange a list of plots via grid.arrange(), but sharing the legend from the first plot
arrange.share.legend <- function(plots, whose.legend = 1, legend.pos = "right",
								 orientation = c("vertical","horizontal"), widths = NULL, heights = NULL, ...) {
	
	nplots <- length(plots)
	
	#print(lapply(plots, class))
	g <- ggplotGrob(plots[[whose.legend]] + theme(legend.position = legend.pos))$grobs
	legend <- g[[ which(sapply(g, function(x) x$name) == "guide-box") ]]
	nolegend <- do.call(arrangeGrob, lapply(plots, function(x) {
		if (inherits(x, "ggplot"))
			x + theme(legend.position = "none")
		else
			return(x)
	}))
	
	orientation <- match.arg(orientation)
	if (orientation == "vertical") {
		lheight <- sum(legend$height)
		if (is.null(heights) || length(heights) != nplots)
			heights <- rep(1, nplots)
		heights <- heights/sum(heights)
		arrangeGrob(nolegend, legend, ncol = 1,
					heights = unit.c(unit(1, "npc") - lheight, lheight))
	}
	else if (orientation == "horizontal") {
		lwidth <- sum(legend$width)
		if (is.null(widths) || length(widths) != nplots)
			widths <- rep(1, nplots)
		widths <- widths/sum(widths)
		arrangeGrob(nolegend, legend, nrow = 1, ...,
					widths = unit.c(unit(1, "npc") - lwidth, lwidth))
	}
	
}

## reorder a factor within groups, so final ordering is nice for a plot axis with panels
reorder.by.group <- function(x, g, v, fn, ...) {
	
	ll <- split(x, g)
	vv <- split(v, g)
	
	rez <- list()
	ymax <- 0
	for (i in seq_along(vv)) {
		o <- order(vv[[i]])
		rez[[i]] <- o+ymax
		ymax <- ymax + max(o)
	}
	
	rez <- unlist(rez)
	print(cbind(unsplit(ll, g), rez))
	reorder(factor(unsplit(ll, g)), rez)
	
}