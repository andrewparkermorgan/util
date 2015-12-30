## ms.R
## read output of runs of Hudson's <ms> simulator

library(pegas)
library(ape)

read.ms <- function(ff, ...) {
	
	if (!file.exists(ff))
		stop("Can't find file <", ff, ">.")
	
	ll <- readLines(ff)
	starts <- grep("//", ll)+3
	pops <- vector("list", length(starts))
	for (i in seq_along(starts)) {
		sites <- unlist(strsplit(gsub("positions: ", "", ll[ starts[i]-1 ]), " "))
		end <- which(grepl("^\\s*$", ll) & 1:length(ll) > starts[i])
		samples <- ll[ (starts[i]):(end[1]-1) ]
		pops[[i]] <- t(sapply(strsplit(samples, ""), as.numeric))
		attr(pops[[i]], "sites") <- as.numeric(sites)
	}
	
	class(pops) <- c("ms", class(pops))
	return(pops)
	
}

sites.ms <- function(pops, ...) {
	lapply(pops, function(x) attr(x, "sites"))
}

nsites.ms <- function(pops, ...) {
	sapply(sites.ms(pops), length)
}

.as.nucleotides <- function(pop, ...) {
	
	alphabet <- c("a","g")
	if (!is.null(attr(pop, "sites")))
		sites <- attr(pop, "sites")
	else
		sites <- NULL
	
	alphabet <- c("a","g")
	x <- matrix(alphabet[ pop+1 ], nrow = nrow(pop), ncol = ncol(pop))
	attr(x, "sites") <- sites
	if (inherits(pop, "ms"))
		class(x) <- c("ms", class(x))
	return(x)
	
}

theta.ms <- function(pop, i = TRUE, ...) {

	x <- .as.nucleotides(pop)
	sites <- attr(x, "sites")
	s <- length(sites[i])
	pegas::theta.s(s, nrow(x))
	
}

nucdiv.ms <- function(pop, i, ...) {
	
	x <- .as.nucleotides(pop)
	x <- x[ ,i, drop = FALSE ] 
	#print(x)
	pegas::nuc.div(as.DNAbin(x))
	
}