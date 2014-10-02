read.dotplot <- function(ff, offset = c(0,0), ...) {
	
	df <- read.table(ff)
	cols <- colnames(df)
	if (ncol(df) == 11)
		## assume this was a self-alignment and is missing 'name2' column
		cols <- c("name1","start1","end1","strand1","start2","end2","strand2",
							"ident","ident.pct","continuity","continuity.pct")
	else
		## assume this was an alignment using my standard lastz output style
		cols <- c("name1","start1","end1","strand1","name2","start2","end2","strand2",
							"ident","ident.pct","continuity","continuity.pct")
	colnames(df) <- cols
	
	if (length(offset) < 2 & length(offset))
		offset[2] <- offset[1]
	
	df$ident.pct <- as.numeric(gsub("%","", df$ident.pct))/100
	df$continuity.pct <- as.numeric(gsub("%","", df$continuity.pct))/100
	df$width <- with(df, end1-start1)
	df$start1 <- df$start1 + offset[1]
	df$end1 <- df$end1 + offset[1]
	df$start2 <- df$start2 + offset[2]
	df$end2 <- df$end2 + offset[2]
	
	return(df)
	
}

ggdotplot <- function(dp, xlim = NULL, ylim = NULL, ...) {
	
	require(ggplot2)
	
	if (is.null(ylim) & !is.null(xlim))
		ylim <- xlim
	if (is.null(xlim) & !is.null(ylim))
		xlim <- ylim
		
	ggplot(dp) +
		geom_segment(aes(x = start1, xend = end1, y = start2, yend = end2, alpha = ident.pct)) +
		geom_abline(slope = 1, intercept = 0, col = "grey70", lty = "dashed") +
		scale_x_continuous(label = function(x) paste(x/1e6)) +
		scale_y_continuous(label = function(x) paste(x/1e6)) +
		scale_alpha_continuous("sequence identity") +
		guides(colour = guide_legend(title = NULL)) +
		coord_fixed(xlim = xlim, ylim = ylim) +
		theme_bw() +
		xlab("\nposition (Mbp)") + ylab("position (Mbp)\n")
	
}