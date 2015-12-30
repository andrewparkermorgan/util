#library(MUGAExampleData)
library(lme4)
library(Matrix)
library(DOQTL)
library(plyr)
library(reshape2)

source("~/lib/util/R/mouse.R")

## extract a slice of DOQTL's haplotype-probability array into a dataframe
.extract.locus <- function(probs, locus, ind = TRUE, ...) {
	x <- probs[ ind,,locus, drop = TRUE ]
	x[ is.na(x) | !is.finite(x) | is.nan(x) ] <- 0
	as.data.frame(x)
}

## construct dataframe for use with lme4 from DOQTL's input objects
.make.model.input <- function(formula, pheno, probs, K, locus, ...) {
	
	fm <- lFormula(formula = as.formula(formula), data = pheno,
					control = lmerControl(check.nobs.vs.nlev = "ignore",
										  check.nobs.vs.nRE = "ignore"))
	df <- fm$fr
	keep <- intersect(intersect(rownames(df), rownames(probs)), rownames(K))
	df <- df[ keep,, drop = FALSE ]
	geno <- .extract.locus(probs, locus, keep)
	df <- cbind(df, geno)
	df$animal <- seq_len(nrow(df))
	rez <- list(df = df, k = K[ keep,keep ], formula = formula)
	return(rez)
	
}

## fit null and full models as a locus with lme4, using result from .make.model.input()
.fit.lmm.locus <- function(x, geno = LETTERS[1:8], intercept = TRUE, family = NULL, ...) {
	
	formula <- x$formula
	data <- x$df
	K <- x$k
	
	fm <- as.formula(formula)
	tm <- terms(fm)
	tl <- as.character(fm)
	if (!attr(tm, "response"))
		tl <- c("pheno", tl[-1])
	else
		tl <- tl[-1]
	
	#fixef <- tl[ !grepl("|", tl) ]
	
	if (!intercept)
		formula <- paste(tl[1], "~ 0 +", paste(tl[-1], collapse = " + "))
	else
		formula <- paste(tl[1], "~ ", paste(tl[-1], collapse = " + "))
	
	f1 <- paste0(formula, " + ", paste(geno, collapse = " + "), " + (1|animal)")
	f0 <- paste0(tl[1], " ~ 1 + ", paste(tl[-1], collapse = " + "), " + (1|animal)")
	ff <- list(reduced = f0, full = f1)
	
	.lmer.animal.model <- function(formula) {
		
		if (is.null(family))
			mf1 <- lFormula(formula = as.formula(formula), data = data,
							control = lmerControl(check.nobs.vs.nlev = "ignore",
												  check.nobs.vs.nRE = "ignore"))
		else
			mf1 <- glFormula(formula = as.formula(formula), data = data, family = family,
							control = glmerControl(check.nobs.vs.nlev = "ignore",
												  check.nobs.vs.nRE = "ignore"))
		
		relmat <- list(animal = K)
		relfac <- relmat
		flist <- mf1$reTrms[["flist"]] # list of grouping factors
		Ztlist <- mf1$reTrms[["Ztlist"]] # list of transpose of the sparse model matrices
		stopifnot(all(names(relmat) %in% names(flist)))
		asgn <- attr(flist, "assign")
		for(i in seq_along(relmat)) {
			tn <- which(match(names(relmat)[i], names(flist)) == asgn)
			if(length(tn) > 1)
				stop("a relationship matrix must be associated",
					 " with only one random effects term", call.=FALSE)
			relmat[[i]] <- Matrix(relmat[[i]], sparse=TRUE)
			relfac[[i]] <- chol(relmat[[i]])
			Ztlist[[i]] <- relfac[[i]] %*% Ztlist[[i]]
		}
		mf1$reTrms[["Ztlist"]] <- Ztlist
		mf1$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
		
		if (is.null(family)) {
			devianceFunction <- do.call(mkLmerDevfun, mf1)
			optimizerOutput <- optimizeLmer(devianceFunction)		
		}	
		else {
			devianceFunction <- with(mf1, mkGlmerDevfun(fr, X, reTrms, family = family))
			optimizerOutput <- optimizeGlmer(devianceFunction)
		}
		
		fit.lme4 <- mkMerMod(rho=environment(devianceFunction),
							 opt=optimizerOutput,
							 reTrms=mf1$reTrms,
							 fr=mf1$fr)
		
		return(fit.lme4)
		
	}
	
	rez <- lapply(ff, .lmer.animal.model)
	class(rez) <- c("doqtl.mm", class(rez))
	return(rez)
	
}

## get rough estimate of "marginal" h^2
## for details see Nakagawa 2013 http://dx.doi.org/10.1111/j.2041-210x.2012.00261.x
h2 <- function(x, df = 8, ...) {
	
	x <- x$full
	vc <- c(unlist(VarCorr(x)), attr(VarCorr(x),"sc")^2)
	
	X <- x@pp$X
	asgn <- attr(X, "assign")
	ncol <- max(asgn)
	i <- which(asgn >= (ncol - df + 1))
	yhat <- X[ ,i ] %*% fixef(x)[i]
	sigma.f <- as.vector(var(yhat))
	
	R2 <- (sigma.f + vc["animal"])/(sum(vc)+sigma.f)
	return(unname(R2))
	
}

## extract fixed effects and their standard errors
fixef.se <- function(x, ...) {
	fe <- fixef(x)
	se <- sqrt(diag(vcov(x)))
	data.frame(coef = names(fe), beta = fe, se = se)
}

## force Sanger-like to UCSC-like chromosome names
.fix.chrs <- function(x, ...) {
	
	x <- as.character(x)
	chrs <- paste0("chr", c(1:19,"X","Y","M"))
	
	if (!any(grepl("^chr", x)))
		x <- paste0("chr", x)
	
	factor(x, levels = chrs)
	
}

## convert a 'doqtl' object (result of DOQTL::scanone()) to a list of dataframes,
##	one for LOD scores and one for coefficient estimates
.melt.doqtl <- function(x, ...) {
	
	lod <- ldply(x$lod, rbind)
	colnames(lod)[1:4] <- c("chr.type","marker","chr","pos")
	lod$pos <- floor(1e6*lod$pos)
	rownames(lod) <- as.character(lod$marker)
	lod$chr <- .fix.chrs(lod$chr)
	
	snps <- unique(lod[ ,c("marker","chr","pos","cM") ])
	
	betas <- ldply(x$coef, function(y) melt(y[ ,-1, drop = FALSE ]))
	colnames(betas) <- c("chr.type","marker","coef","beta")
	betas <- merge(betas, snps)
	
	rez <- list(betas = arrange(betas, chr.type, chr, pos, coef),
				lod = arrange(lod, chr.type, chr, pos))
	
	return(rez)
	
}

## plot LOD scores along concatenated chromosomes (like R/qtl::plot.scanone())
lodplot <- function(x, show = c("lod","neg.log10.p"), max.lod = Inf, ...) {
	
	if (!inherits(x, "doqtl"))
		stop("Please supply an object of class 'doqtl', the result from DOQTL::scanone().")
	
	show <- match.arg(show)
	
	x <- .melt.doqtl(x)
	o <- x$lod[ ,show ] > max.lod
	no <- sum(o, na.rm = TRUE)
	if (no > 0)
		message("Truncating '", show, "' to ", round(max.lod, 2), "; ",
				no, " markers omitted, in addition to those with NA/NaN values.")
	
	x$lod[ ,show ] <- pmin(x$lod[ ,show ], max.lod)
	p <- ggmanhattan(x$lod) +
		geom_line(aes_string(y = show, group = "chr"))
	
	if (show == "lod")
		p <- p + ylab("LOD score\n")
	else if (show == "neg.log10.p")
		p <- p + ylab("-log10( p-value )\n")
	
	return(p)
	
}


## plot allele effects (= model coefficients) and LOD scores in parallel along single chromosome
coefplot <- function(x, chroms = "chr1", rows, max.beta = Inf, coefs = LETTERS[1:8], ...) {
	
	if (!inherits(x, "doqtl"))
		stop("Please supply an object of class 'doqtl', the result from DOQTL::scanone().")
	
	colnames(x$coef$A)[1] <- "A"
	nm <- colnames(x$coef$A)
	#if (nchar(nm[2]) == 2) {
	#	nm <-substr(nm, 1, 1)
	#}
	#colnames(x$coef$A)[1] <- nm
	
	x <- .melt.doqtl(x)
	o <- x$betas$beta > max.beta
	no <- sum(o, na.rm = TRUE)
	x$betas$beta <- pmin(x$betas$beta, max.beta)
	if (no > 0)
		message("Truncating '", show, "' to ", round(max.beta, 2), "; ",
				no, " markers omitted, in addition to those with NA/NaN values.")
	
	toplot <- subset(x$betas, chr == chroms & coef %in% coefs)[ ,c("chr","marker","pos","coef","beta") ]
	toplot$panel <- "coefficients"
	toplot$lod <- NA
	
	lods <- x$lod[ unique(toplot$marker),c("chr","marker","pos","lod"), drop = FALSE ]
	lods$beta <- NA
	lods$coef <- NA
	lods$panel <- "LOD scores"
	
	toplot <- rbind(toplot, lods)
	if (!missing(rows)) {
		i <- eval(substitute(rows), toplot, parent.frame())
		toplot <- toplot[i,]
	}
	
	ggplot(toplot) +
		geom_line(aes(x = pos, y = beta, colour = coef)) +
		geom_line(aes(x = pos, y = lod)) +
		scale_x_genome() +
		scale_colour_manual("founder", values = CC.COLORS[1:8], labels = setNames(cc.strains,LETTERS[1:8])) +
		facet_grid(panel ~ ., scale = "free_y") +
		xlab(paste("\n", chroms, "position (Mbp)")) + ylab("")
	
}

## plot estimated allele effects at a single locus, using output from .fit.mm.locus()
locusplot <- function(x, geno = LETTERS[1:8], ...) {
	
	betas <- fixef.se(x$full)
	ggplot(subset(betas, coef %in% geno)) +
		geom_pointrange(aes(x = coef, colour = coef,
							y = beta, ymin = beta-1.96*se, ymax = beta+1.96*se)) +
		scale_colour_manual("founder", values = CC.COLORS[1:8], labels = setNames(cc.strains,LETTERS[1:8]))
	
}

## plot haplotype reconstruction from the 36-state probability matrix for a single individual
plot.doqtl.genome <- function(probs, snps, colours = NULL, ...) {
	
	## get max-likelilhood states
	message("Determining maximum-likelihood diplotype states...")
	states <- colnames(probs)
	maxstates <- apply(probs, 1, which.max)
	haps <- ldply(strsplit(states[ maxstates ], ""))
	colnames(haps)[1:2] <- c("hap1","hap2")
	haps$marker <- rownames(probs)
	haps <- merge(haps, snps, all.x = TRUE)
	
	## do pseudo-phasing
	pseudophase <- function(d) {
		
		message("pseudo-phasing chromosome ", d$chr[1], " ...")
		if (nrow(d) > 1) {
			for (i in 2:nrow(d)) {
				if (d$hap1[i] == d$hap2[i-1]) {
					#message("\tswapping haplotypes at marker ", i)
					x <- d$hap2[i]
					d$hap2[i] <- d$hap1[i]
					d$hap1[i] <- x
				}
				else if (d$hap2[i] == d$hap1[i-1]) {
					#message("\tswapping haplotypes at marker ", i)
					x <- d$hap1[i]
					d$hap1[i] <- d$hap2[i]
					d$hap2[i] <- x
				}
			}
			return(d)
		}
		else {
			return(d)
		}
		
	}
	
	rle.to.blocks <- function(run) {
		j <- 0
		starts <- integer(0)
		ends <- integer(0)
		values <- character(0)
		for (i in seq_along(run$lengths)) {
			starts <- c(starts, j+1)
			j <- j+run$lengths[i]
			ends <- c(ends, j)
			values <- c(values, run$values[i])
		}
		data.frame(start = starts, end = ends, value = values)
	}
	
	haps <- arrange(haps, chr, pos)
	phased <- ddply(haps, .(chr), pseudophase)
	
	phased.m <- ddply(phased, .(chr), function(d) {
		diplo <- paste0(d$hap1, d$hap2)
		blocks <- rle.to.blocks(rle(diplo))
		blocks$start <- d[ blocks$start,"pos" ]
		blocks$end <- d[ blocks$end,"pos" ]
		diplo <- ldply(strsplit(as.character(blocks$value), ""))
		blocks$hap1 <- diplo[,1]
		blocks$hap2 <- diplo[,2]
		blocks$value <- NULL
		melt(blocks, measure.vars = c("hap1","hap2"))
	})
	
	## make plot
	message("generating plot...")
	phased.m$updown <- c(-1,1)[ as.numeric(phased.m$variable) ]*0.5
	chrs <- unique(as.numeric(phased.m$chr))
	chr.labs <- gsub("^chr", "", levels(phased.m$chr)[chrs])
	p <- ggplot(phased.m) +
		geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(chr)-updown*0.05, ymax = as.numeric(chr)-updown*0.8, fill = value)) +
		scale_y_continuous(breaks = chrs, labels = chr.labs, trans = "reverse") +
		scale_x_continuous("\nposition (Mbp)", labels = function(x) x/1e6) +
		scale_fill_discrete("strain")
	
	if (!is.null(colours))
		p <- p + scale_fill_manual("strain", values = colours)
	
	return(p)
	
}

