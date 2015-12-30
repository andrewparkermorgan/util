## plot phylogeny annotated with (known) taxa assignments

pretty.phylo <- function(tree, bywhat, colors = NULL, widths = c(1,3), pch = 19, cex = 1,
						 label = NULL, show.tip.label = FALSE, type = "fan", fill = NULL,
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
	colors.final <- colors[ as.character(bywhat) ]
	fill.final <- colors.final
	if (!is.null(fill)) {
		fill.final <- fill[ as.character(bywhat) ]
	}
	colors.final[ is.na(fill.final) | is.na(colors.final) | is.na(pch) ] <- "grey50"
	fill.final[ is.na(fill.final) | is.na(colors.final) | is.na(pch) ] <- "grey50"
	pch[ is.na(pch) ] <- 19
	
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
	plot(tree, type = type, show.tip.label = show.tip.label, cex = cex, ...)
	tiplabels(pch = pch, ..., col = colors.final, bg = fill.final)
	
}

plot.with.bootstrap <- function(tree, cex = 1, ...) {
	
	plot(tree, cex = cex, ...)
	nodelabels(tree$node.label, bg = NA, frame = "none", cex = 0.75*cex, col = "grey50", adj = 0)
	
}

plot.omega <- function(aln, ...) {
	
	require(reshape2)
	require(ggplot2)
	
	.kaks <- calc.omega(aln)
	ks.df <- melt(.kaks)
	colnames(ks.df) <- c("seq1","seq2","kaks")
	ks.df$kaks[ !is.finite(ks.df$kaks) | abs(ks.df$kaks) > 10 | ks.df$seq1 == ks.df$seq2 ] <- NA
	ks.df$seq1 <- factor(ks.df$seq1, levels = colnames(.kaks))
	ks.df$seq2 <- factor(ks.df$seq2, levels = colnames(.kaks))
	ggplot(ks.df, aes(x = seq1, y = seq2, fill = kaks)) +
		geom_tile() +
		scale_fill_gradient("Ka/Ks") +
		coord_equal() +
		theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
						   axis.title = element_blank())
	
}

## draw phylogenetic trees with ggplot2
## 	with thanks to the phyloseq package (https://github.com/joey711/phyloseq)
ggplot.phylo <- function(phy, data = NULL, type = "default", label.tips = TRUE, label.nodes = FALSE, label.by = "taxa",
						 use.edge.length = TRUE, spread = 0.1, lab.offset = 0.02, lab.size = 3, justify = FALSE,
						 root.age = NULL, expand = c(0.2, 0), ...) {
	
	require(ggplot2)
	require(scales)
	require(grid)
	require(ape)
	require(plyr)
	
	## function to return node IDs for tips, *in current tree order*
	.get.tips <- function(t, ...) {
		t$edge[ !(t$edge[,2] %in% t$edge[,1]),2 ]
	}
	
	## define some function calls to ape internals
	ape_node_height <- function(Ntip, Nnode, edge, Nedge, yy){
		.C(ape:::node_height, PACKAGE = "ape",
		   as.integer(Ntip), as.integer(Nnode),
		   as.integer(edge[ ,1 ]), as.integer(edge[ ,2 ]),
		   as.integer(Nedge), as.double(yy))[[6]]
	}
	ape_node_depth <- function(Ntip, Nnode, edge, Nedge, node.depth) {
		.C(ape:::node_depth, PACKAGE = "ape",
		   as.integer(Ntip), as.integer(Nnode),
		   as.integer(edge[ ,1 ]), as.integer(edge[, 2]),
		   as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[6]]
	}
	ape_node_depth_edge_length <-function(Ntip, Nnode, edge, Nedge, edge.length) 
		.C(ape:::node_depth_edgelength, PACKAGE = "ape",
		   as.integer(Ntip), as.integer(Nnode),
		   as.integer(edge[, 1]), as.integer(edge[,2]),
		   as.integer(Nedge), as.double(edge.length), double(Ntip + Nnode))[[7]]
	
	## order the tree
	z <- reorder.phylo(phy, order = "postorder")
	
	## prepare some stuff for calculating edge positions
	Nedge <- nrow(phy$edge)[1]
	Nnode <- phy$Nnode
	Ntip <- length(phy$tip.label)
	ROOT <- Ntip + 1
	TIPS <- phy$edge[(phy$edge[, 2] <= Ntip), 2]
	NODES <- (ROOT):(Ntip + Nnode)
	nodelabels <- phy$node.label
	
	## use ape internals to compute start and end positions of edges
	if (use.edge.length)
		xx <- ape_node_depth_edge_length(Ntip, Nnode, z$edge, Nedge, z$edge.length)
	else {
		z$edge.length <- rep(1, length(z$edge.length))
		#xx <- ape_node_depth(Ntip, Nnode, z$edge, Nedge, 1) - 1
		#xx <- max(xx) - xx
		xx <- ape_node_depth_edge_length(Ntip, Nnode, z$edge, Nedge, z$edge.length)
	}
	yy <- numeric(Ntip + Nnode)
	yy[TIPS] <- 1:Ntip
	yy <- ape_node_height(Ntip, Nnode, z$edge, Nedge, yy)
	eps <- max(xx)*lab.offset
	
	tdf <- data.frame(phy$edge, taxa = NA)
	tl <- phy$tip.label[ .get.tips(phy) ]
	tdf$taxa[ tdf$X2 <= Ntip ] <- tl
	tdf$xleft <- xx[ tdf$X1 ]
	tdf$xright <- xx[ tdf$X2 ]
	tdf$y <- yy[ tdf$X2 ]
	tdf$eps <- eps
		
	edf <- ddply(tdf, .(X1), summarize,
				 x = xleft[1], vmin = min(y), vmax = max(y))
	
	## position internal nodes
	if (!is.null(phy$node.label)) {
		tdf$x <- NA
		tdf$x[ tdf$X2 > ROOT ] <- tdf$xright[ tdf$X2 > ROOT ]
		tdf$label <- NULL
		tdf$label[ tdf$X2 > ROOT ] <- as.numeric(phy$node.label[-1])
	}
	
	## set angles for text
	ymax <- (1+spread)*max(tdf$y, na.rm = TRUE)
	tdf$hj <- 0
	tdf$angle <- with(tdf, -90-360/ymax * y)
	idx <- tdf$angle > -270 & tdf$angle < -90
	tdf$angle[idx] <- 180+tdf$angle[idx]
	tdf$hj[!idx] <- 1
	
	## merge in node data
	if (!is.null(data))
		tdf <- merge(tdf, data, all.x = TRUE)
	nodes <- subset(tdf, !is.na(taxa))
	
	if (type == "fan") {
		if (!justify)
			node.aes <- aes_string(angle = "angle", hjust = "hj", label = label.by, x = "xright + eps")
		else
			node.aes <- aes_string(angle = "angle", hjust = "hj", label = label.by, x = "max(xright)+eps")
	}
	else {
		if (justify)
			xpos = "max(xright)+eps"
		else
			xpos = "xright + eps"
		node.aes <- aes_string(label = label.by, hjust = 0, x = xpos)
	}
		
	## make base plot
	#print(tdf)
	#print(nodes)
	p <- ggplot(tdf, aes(x = xright, y = y)) +
		geom_segment(aes(x = xleft, xend = xright, y = y, yend = y)) +
		geom_segment(data = edf, aes(x = x, xend = x, y = vmin, yend = vmax)) +
		scale_y_continuous(limits = c(0, ymax), expand = c(0, 1))
	
	if (label.tips) {
		p <- p + geom_text(data = tdf[ !is.na(tdf[,label.by]), ], node.aes, size = lab.size)
		if (justify)
			p <- p + geom_segment(data = nodes,
								  aes(x = xright, xend = max(xright)-0.01*max(xright), y = y, yend = y),
								  colour = "grey70", size = 0.1)
	}
	
	if (label.nodes) {
		node.aes <- aes_string(x = "xright + 1.8*eps", label = "label", size = lab.size, hjust = 0)
		p <- p + geom_text(data = tdf[ !is.na(tdf[,"label"]), ], node.aes, size = 0.75*lab.size, col = "grey50")
	}
	
	if (!is.null(root.age))
		p <- p + scale_x_continuous(breaks = seq(0, max(tdf$xright, na.rm = TRUE), length.out = 5),
									labels = function(x) root.age - root.age*x/max(tdf$xright, na.rm = TRUE),
									expand = expand)
	
	if (type == "fan")
		p <- p + coord_polar("y")
	
	return(p)
	
}

## a useful default theme for ggplot phylogeny
theme_phylo <- function(...) {

	thm <- theme_classic(...)
	thm + theme(axis.ticks = element_blank(), axis.text = element_blank(),
				axis.title = element_blank(), axis.line = element_blank())
	
}

## like above, but without margins: synonymous with Hadley's theme_nothing()
theme_phylo_nolegend <- function(base_size = 12, base_family = "Helvetica") {
	
		theme_bw(base_size = base_size, base_family = base_family) %+replace%
			theme(
				rect             = element_blank(),
				line             = element_blank(),
				text             = element_blank(),
				axis.ticks.margin = unit(0, "lines")
			)
		
}