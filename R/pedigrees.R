## check pedigrees

ped.to.graph <- function(x, effect = NULL, ...) {
	## from synbreed::plot.pedigree()
	require(igraph)
	
	if (any(class(x) == "gpData")) 
		pedigree <- x$pedigree
	else pedigree <- x
	if (min(pedigree$gener) != 0) 
		pedigree$gener <- pedigree$gener - min(pedigree$gener)
	if (!is.null(effect) & length(effect) != nrow(pedigree)) 
		stop("length of effect does not equal nrow(pedigree)")
	pedigree[!(pedigree$Par1 %in% pedigree$ID), ]$Par1 <- 0
	pedigree[!(pedigree$Par2 %in% pedigree$ID), ]$Par2 <- 0
	relations <- rbind(as.matrix(pedigree[pedigree$Par1 != 0, 
										  c("Par1", "ID")], ncol = 2), as.matrix(pedigree[pedigree$Par2 != 
										  													0, c("Par2", "ID")], ncol = 2))
	ped.graph <- graph.data.frame(relations, directed = TRUE, 
								  pedigree)
}

layout.pedigree <- function(x, effect = NULL, ...) {
	
	if (any(class(x) == "gpData")) 
		pedigree <- x$pedigree
	else pedigree <- x
	if (min(pedigree$gener) != 0) 
		pedigree$gener <- pedigree$gener - min(pedigree$gener)
	if (!is.null(effect) & length(effect) != nrow(pedigree)) 
		stop("length of effect does not equal nrow(pedigree)")
	
	gener <- pedigree$gener
	n <- nrow(pedigree)
	pos <- matrix(data = NA, nrow = n, ncol = 2)
	pos[, 2] <- max(gener) - gener
	if (is.null(effect)) 
		pos[, 1] <- order(gener, partial = order(pedigree$ID, 
												 decreasing = TRUE)) - cumsum(c(0, table(gener)))[gener + 
												 												 	1]
	else pos[, 1] <- effect
	myscale <- function(x) {
		if (length(x) == 1) 
			x <- 0
		else {
			x <- unlist(scale(x))
		}
		return(x)
	}
	if (is.null(effect)) 
		pos[n:1, 1] <- unlist(tapply(pos[, 1], pos[, 2], myscale))
	
	return(pos)
	
}
