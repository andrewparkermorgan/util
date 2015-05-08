## Convenience functions for constructing graphs or exploring graph structure.
## All these are thin wrappers of igraph functions.

library(igraph)

## get integer vertex ID from vertex name
vid <- function(vname, g, ...) {
	match(vname, V(g)$name)
}

## does there exist any path between node pairs <v1>,<v2>?
exists.path <- function(g, v1, v2, ...) {
	
	if (length(v1) != length(v2))
		stop("")
	
	plen <- shortest.paths(g, v1, v2, ...)
	is.finite(plen[ cbind(v1, v2) ])
	
}

## enumerate nodes which can be reached from a focal node
## returns list (so can specify >1 focal node)
connections <- function(g, focal, ...) {
	
	sp <- shortest.paths(g)
	
	.get.connections <- function(n) {
		V(g)$name[ is.finite(sp[ n, ]) ]
	}
	
	lapply(focal, .get.connections)
	
}

## get largest conencted component of a graph
biggest.cluster <- function(g, ...) {
	
	cl <- clusters(g)
	i <- which.max(cl$csize)
	induced.subgraph(g, cl$membership == i)
	
}

## construct weighted graph from tuples of (node1, node2, weight)
graph.from.edges <- function(v1, v2, weights = 1, weighted = TRUE, ...) {
	
	require(plyr)
	
	v1 <- as.character(v1)
	v2 <- as.character(v2)
	vv <- union(v1, v2)
	A <- matrix(0, nrow = length(vv), ncol = length(vv))
	rownames(A) <- colnames(A) <- vv
	
	df <- data.frame(v1 = v1, v2 = v2, weights = weights)
	ddf <- ddply(df, .(v1, v2), summarize,
				 w = sum(weights), n = length(weights))
	
	if (weighted)
		ecol <- "w"
	else
		ecol <- "n"
	
	A[ with(ddf, cbind(as.character(v1), as.character(v2))) ] <- ddf[ ,ecol ]
	
	if (!weighted)
		weighted <- NULL
	
	graph.adjacency(A, weighted = weighted)
	
}