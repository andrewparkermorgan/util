library(ape)
library(phytools)

theme_phylo_blank <- function(family = "Helvetica", size = 12, face = "italic", colour = "black")
{
	theme_bw() %+replace%
		theme(
			text = element_text(family = family, size = size, face = face, hjust = 0, vjust = 0, angle = 0, lineheight = 1, colour = colour),
			panel.border = element_blank(),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			plot.background = element_blank(),
			strip.background = element_blank(),
			axis.text = element_blank(),
			axis.ticks = element_blank(),
			axis.title = element_blank(),
			axis.ticks.margin = unit(0, "lines")
		)
}

get.tips <- function(t, ...) {
	t$edge[ !(t$edge[,2] %in% t$edge[,1]),2]
}

find.tip <- function(t, name, ...) {
	
	i <- t$tip.label == name
	if (sum(i))
		return( which(i) )
	else
		return(NA)
	
}

grep.tips <- function(t, pattern, ...) {
	
	grep(pattern, t$tip.label, ...)
	
}

get.parent <- function(t, node, ...) {
	
	if (is.null(t$tips))
		t$tips <- get.tips(t)
	
	return( t$edge[ t$edge[,2] == node,1 ] )
	
}

is.root <- function(t, node, ...) {
	!(node %in% t$edge[,2])
}

traverse.up <- function(start, tree, time = Inf, dist = NULL, verbose = FALSE, ...) {
	
	t <- tree
	if (!is.rooted(t))
		stop("Can't traverse backward in time on an unrooted tree.")
	
	if (is.null(dist))
		dist <- dist.nodes(t)
		
	if (!is.root(t, start)) {
		time <- abs(time)
		last <- 0
		elapsed <- 0
		p <- start
		prev <- start
		while (elapsed < time & !is.root(t, start)) {
			p <- get.parent(t, start, time - elapsed)
			step <- dist[p,start]
			last <- elapsed
			elapsed <- elapsed + step
			prev <- start
			start <- p
			if (verbose) {
				cat("Starting at [", start, "]", "with", time, "distance left\n")
				cat("distance to", p, "was", step, "; we've spent", elapsed, "\n")
			}
		}
		return( matrix(c(prev, p, last, elapsed), nrow = 2) )
	}
	else {
		return(NA)
	}
	
}