library(ggplot2)
library(igraph)

## use ggplot to make a nice-looking plot of an igraph object
## <aug> = metadata for nodes, <merge.by> = field in <aug> on which to join
ggplot.igraph <- function(g, layout = layout.kamada.kawai,
						  aug = NULL, merge.by = NULL, ...) {
	
	## use canned igraph layout function to get xy-coords for nodes
	nodes <- as.data.frame(layout(g))
	colnames(nodes) <- c("x","y")
	
	## add 'vertex attributes' (node metadata)
	for (a in list.vertex.attributes(g)) {
		nodes[ ,a ] <- get.vertex.attribute(g, a)
	}
	
	## add external node metadata, if provided
	if (!is.null(aug))
		if (!is.null(merge.by))
			nodes <- merge(nodes, aug, by.x = "name", by.y = merge.by, all.x = TRUE)
	else
		nodes <- merge(nodes, aug, all.x = TRUE)
	
	## use node names as rownames
	rownames(nodes) <- nodes$name
	print(head(nodes))
	
	## make separate dataframe for edges, respecting weights if present
	el <- get.edgelist(g)
	edges <- data.frame(from.x = nodes[ el[,1],"x" ], from.y = nodes[ el[,1],"y" ],
						to.x = nodes[ el[,2],"x" ], to.y = nodes[ el[,2],"y" ])
	if (is.weighted(g)) {
		edges$weight <- E(g)$weight
		edge.aes <- aes(x = from.x, xend = to.x, y = from.y, yend = to.y, alpha = weight)
	}
	else {
		edge.aes <- aes(x = from.x, xend = to.x, y = from.y, yend = to.y)
	}
	print(head(edges))
	
	## return skeleton plot; user can customize by adding further layers
	ggplot(nodes, aes(x = x, y = y), ...) +
		geom_segment(data = edges, edge.aes)
	
}

## blank canvas for a graph layout, but preserving space for legend(s)
theme_igraph <- function(...) {
	
	theme_classic(...) +
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
			  axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
			  axis.title = element_blank())
	
}