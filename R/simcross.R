## auxiliary functions for closed-population breeding simulations
## largely adapted from kwbroman's R/simcross package <https://github.com/kbroman/simcross/>

library(simcross)
library(plyr)

## simulate an outbred pedigree from a given number of founders
## NB: founders are automatically expanded so that one of each sex occurs in 0th generation
outbred.from.founders <- function(nfounders = 8, npairs = nfounders/2, ngen = 2,
								  design = c("nosib","random"), ...) {
	
	stopifnot(nfounders %% 2 == 0)
	stopifnot(npairs > 0)
	design <- match.arg(design)
	
	## initial generation
	id <- c(1:nfounders, nfounders+(1:nfounders))
	mom <- dad <- rep(0, nfounders*2)
	sex <- rep(c(0,1), nfounders)
	gen <- rep(0, nfounders*2)
	
	dads <- id[ sample(which(sex == 0), npairs, replace = (npairs > nfounders)) ]
	moms <- id[ sample(which(sex == 1), npairs, replace = (npairs > nfounders)) ]
	
	for (i in 1:ngen) {
		
		while (design == "nosib" && any(dads - moms == 1)) {
			# sample until no sibs
			dads <- sample(dads)
		}
		kids <- 1:(npairs*2)+max(id)
		id <- c(id, kids)
		mom <- c(mom, rep(moms, each = 2))
		dad <- c(dad, rep(dads, each = 2))
		sex <- c(sex, rep(c(0,1), npairs))
		gen <- c(gen, rep(i, npairs*2))
		
		moms <- kids[ seq(1, length(kids), by = 2) ]
		dads <- kids[ seq(2, length(kids), by = 2) ]
		
	}
	
	cbind(id = id, mom = mom, dad = dad, sex = sex, gen = gen)
	
}

## add one more outbreeding generation to an existing pedigree
draw.next.generation <- function(ped, npairs = 20, design = c("nosib","random"), ...) {
	
	design <- match.arg(design)
	
	last.id <- max(ped[,1])
	last.gen <- max(ped[,5])
	
	ids <- 1:(2*npairs)+last.id
	moms <- ped[ ped[,5] == last.gen & ped[,4] == 0, 1 ]
	dads <- ped[ ped[,5] == last.gen & ped[,4] == 1, 1 ]
	moms <- sample(moms, npairs, replace = (npairs > length(moms)))
	dads <- sample(dads, npairs, replace = (npairs > length(dads)))
	
	while (design == "nosib" && any(dads - moms == 1)) {
		# sample until no sibs
		dads <- sample(dads)
	}
	
	mom <- rep(moms, each = 2)
	dad <- rep(dads, each = 2)
	sex <- rep(c(0,1), npairs)
	gen <- rep(last.gen+1, 2*npairs)
	
	cbind(id = ids, mom = mom, dad = dad, sex = sex, gen = gen)
	
}


## simulate chromosomes for next pedigree chunk, given result up to previous chunk
## pedigree = existing pedigree
## ped.new = next pedigree chunk
## geno = existing chromosomes
sim.next.generation <- function(pedigree, ped.new = NULL, geno = NULL, L = 100, xchr = FALSE,
								m = 10, p = 0, obligate_chiasma = FALSE, ...) {
	
	if (length(L) > 1) { # multiple chromosomes
		
		result <- vector("list", length(L))
		if (is.null(names(L)))
			names(L) <- seq(along=L)
		names(result) <- names(L)
		
		if (is.character(xchr)) # xchr is chromosome names
			xchr <- names(L) %in% xchr
		
		if (is.null(xchr))
			xchr <- rep(FALSE, length(L))
		
		if (length(xchr) == 1) # if single value, apply to all chromosomes
			xchr <- rep(xchr, length(L))
		
		if (length(xchr) != length(L))
			stop("length(xchr) != length(L)")
		
		for (i in seq(along=L))
			result[[i]] <- sim_from_pedigree(pedigree, L[i], xchr[i],
											 m, p, obligate_chiasma)
		return(result)
	}
	
	if (length(unique(pedigree[,1])) != nrow(pedigree))
		stop("IDs must be unique")
	rownames(pedigree) <- pedigree[,1]
	
	start <- 1
	if (is.null(geno) & is.null(ped.new))
		result <- vector("list", nrow(pedigree))
	else {
		start <- nrow(pedigree)+1
		pedigree <- rbind(pedigree, ped.new)
		result <- c(geno, vector("list", nrow(ped.new)))
	}
	
	names(result) <- as.character(pedigree[,1])
	
	if (obligate_chiasma)
		Lstar <- calc_Lstar(L, m, p)
	else
		Lstar <- L
	
	for (i in start:nrow(pedigree)) {
		if(pedigree[i,2]==0 || pedigree[i,3]==0) # founder
			result[[i]] <- create_parent(L, allele=pedigree[i,1])
		else {
			mom <- which(pedigree[,1]==pedigree[i,2])
			dad <- which(pedigree[,1]==pedigree[i,3])
			
			if(length(mom) < 1 || length(dad) < 1) {
				print(pedigree[i,,drop=FALSE])
				stop("parents not found")
			}
			if(mom >= i || dad >= i) {
				print(pedigree[i,,drop=FALSE])
				stop("Pedigree problem: parents follow individual")
			}
			
			result[[i]] <- cross(result[[mom]], result[[dad]],
								 m = m, p = p,
								 xchr = xchr, male = (pedigree[i,4] == 1),
								 obligate_chiasma = obligate_chiasma,
								 Lstar = Lstar)
		}
	}
	
	## return only the chromosomes for the new chunk of the pedigree
	return(result[ start:length(result) ])

}

## starting from an initial pedigree, simulate outbreeding until fixation or loss of specified founder alleles
sim.until.fix <- function(ped, target = 1.0, maxgen = Inf, ...) {
	
	fixed <- FALSE
	lost <- FALSE
	geno <- sim_from_pedigree(ped)
	gen <- max(ped[,5])
	while (!(fixed || lost) && gen < maxgen) {
		ped.next <- draw.next.generation(ped, ...)
		geno.next <- sim.next.generation(ped, ped.next, geno, ...)
		maf <- sum(allele.copies(geno.next, ...))
		boundary <- 2*nrow(ped.next)*target
		fixed <- maf >= boundary
		lost <- maf == 0
		gen <- max(ped.next[,5])
		geno <- c(geno, geno.next)
		ped <- rbind(ped, ped.next)
	}
	
	return( list(ped = ped, geno = geno, gen = gen, fixed = fixed, lost = lost) )
	
}

## get count of copies of specified founder alleles per individual
allele.copies <- function(geno, pos = 50, founders = 1, sexless = TRUE, ...) {
	
	tbl <- get_geno(geno, pos)
	tbl[,1] <- as.numeric(tbl[,1] %in% founders)
	tbl[,2] <- as.numeric(tbl[,2] %in% founders)
	
	rowSums(tbl)
	
}

## summarize allele counts by generation at specified locus, given output of sim.until.fix()
## run = list with elements 'ped'=pedigree, 'geno'=chromosomes (as simulated by simcross::sim_from_pedigree())
summarise.run <- function(run, founders = 1, pos = 50, ...) {
	
	run$ped <- transform(as.data.frame(run$ped), x = allele.copies(run$geno, pos = pos, founders = founders))
	summ <- ddply(run$ped, .(gen), summarize, maf = sum(x)/(2*length(x)) )
	summ$fixed <- run$fixed
	summ$lost <- run$lost
	summ$ngen <- run$gen
	
	return(summ)
	
}