## plot phylogeny annotated with (known) taxa assignments

pretty.phylo <- function(tree, bywhat, colors = NULL, widths = c(1,3), pch = 19, cex = 1,
						 label = NULL, show.tip.label = FALSE,
						 legend.title = NULL, legend.pos = "left", legend.cex = 1, ...) {
	
	require(RColorBrewer)
	
	if(length(bywhat) != Ntip(tree))
		stop(paste("Labelling factor must be the right length:",
				   "tree has", Ntip(tree),"tips; only", length(by), "labels present"))
	
	if (!is.factor(bywhat))
		bywhat <- factor(bywhat)
	else
		bywhat <- droplevels(bywhat)
	
	if (is.null(colors)) {
		colors <- colorRampPalette(brewer.pal(min(9, nlevels(bywhat)), "Set1"))(nlevels(bywhat))
		names(colors) <- levels(bywhat)
	}
	colors <- colors[ names(colors) %in% levels(bywhat) ]
	
	if (is.null(label) | length(label) != Ntip(tree))
		label <- bywhat
	tree$tip.label <- as.character(label)
	
	if (is.null(legend.title))
		legend.title <- deparse(substitute(bywhat))
	
	## set up two-panel layout
	layout(matrix(c(1,2), ncol = 2), widths = widths)
	
	## plot legend
	par(cex = 0.9)
	par(mar=c(1,1,1,1)) 
	plot.new()
	legend(legend.pos, names(colors), col = colors, pch = 19,
		   bty = "n", cex = legend.cex, xjust = 1,
		   title = legend.title, title.adj = 0.15)
	
	## plot tree
	par(mar=c(1,1,1,1))
	plot(tree, type = "fan", show.tip.label = show.tip.label, cex = cex, ...)
	tiplabels(pch = pch, ..., col = colors[ as.character(bywhat) ])
	
}

plot.with.bootstrap <- function(tree, cex = 1, ...) {
	
	plot(tree, cex = cex, ...)
	nodelabels(tree$node.label, bg = NA, frame = "none", cex = 0.75*cex, col = "grey50", adj = 0)
	
}