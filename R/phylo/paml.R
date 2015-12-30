library(ape)
library(seqinr)

.write.paml.input <- function(x, file, nbcol = 6, colsep = " ", colw = 10, indent = NULL, blocksep = 1, ...) {
	
	if (inherits(x, "DNAbin")) 
		x <- as.character(x)
	
	aligned <- TRUE
	if (is.matrix(x)) {
		N <- dim(x)
		S <- N[2]
		N <- N[1]
		xx <- vector("list", N)
		for (i in 1:N) xx[[i]] <- x[i, ]
		names(xx) <- rownames(x)
		x <- xx
		rm(xx)
	}
	else {
		N <- length(x)
		S <- unique(unlist(lapply(x, length)))
		if (length(S) > 1) 
			aligned <- FALSE
	}
	
	if (is.null(names(x))) 
		names(x) <- as.character(1:N)
	if (is.null(indent)) 
		indent <- 10
	if (is.numeric(indent)) 
		indent <- paste(rep(" ", indent), collapse = "")
	
	blocksep <- paste(rep("\n", blocksep), collapse = "")
	
	zz <- file(file, "w")
	on.exit(close(zz))
	if (!aligned) 
		stop("sequences must have the same length for\n interleaved or sequential format.")
	cat(N, " ", S, "I\n", sep = "", file = zz)
	if (nbcol < 0) {
		nb.block <- 1
		nbcol <- totalcol <- ceiling(S/colw)
	}
	else {
		nb.block <- ceiling(S/(colw * nbcol))
		totalcol <- ceiling(S/colw)
	}
	SEQ <- matrix("", N, totalcol)
	for (i in 1:N) {
		X <- paste(x[[i]], collapse = "")
		for (j in 1:totalcol) SEQ[i, j] <- substr(X, 1 + (j - 1) * colw, colw + (j - 1) * colw)
	}
	max.nc <- max(nchar(names(x)))
	fmt <- paste("%-", max.nc + 2, "s", sep = "")
	names(x) <- sprintf(fmt, names(x))
	
	colsel <- if (nb.block == 1) 1:totalcol else 1:nbcol
	for (i in 1:N) {
		cat(names(x)[i], file = zz)
		cat(SEQ[i, colsel], sep = colsep, file = zz)
		cat("\n", file = zz)
	}
	if (nb.block > 1) {
		for (k in 2:nb.block) {
			cat(blocksep, file = zz)
			endcolsel <- if (k == nb.block) totalcol else nbcol + 
				(k - 1) * nbcol
			for (i in 1:N) {
				cat(indent, file = zz)
				cat(SEQ[i, (1 + (k - 1) * nbcol):endcolsel], 
					sep = colsep, file = zz)
				cat("\n", file = zz)
			}
		}
	}
	
}

## test for positive selection on coding DNA with PAML
paml.branch.site.test <- function(aln, tree, ...) {
	
	## branch-site test for positive selection:
	##	H0: model = 2; ns.sites = 2; omega = 1(?); fix.omega = 1
	##	Ha: model = 2; ns.sites = 2; omega = 1(?); fix.omega = 0
	message("Null hypothesis: equal rates foreground vs background branches")
	r0 <- paml.codeml(aln, tree, fix.omega = 1, ..., verbose = FALSE)
	message("Alternative hypothesis: accelerated rates foreground vs background branches")
	r1 <- paml.codeml(aln, tree, fix.omega = 0, ..., verbose = FALSE)
	
	chisq <- abs(r1$lnL - r0$lnL)
	df <- abs(r1$df - r0$df)
	rez <- list(method = "Branch-site test by LRT", data.name = "<PAML::codeml output>",
				statistic = c("-2deltalnL" = chisq), parameter = c("df" = df),
				p.value = pchisq(chisq, df, lower.tail = FALSE))
	class(rez) <- "htest"
	return(rez)
	
}

## wrapper for PAML::codeml; captures parts of result necessary to do model tests
paml.codeml <- function(aln = NULL, tree = NULL, outfile = NULL, verbose = TRUE,
						seqtype = 1, clock = 0, model = 2, ns.sites = 2, omega = 0, fix.omega = 0, ...) {
	
	if (is.null(outfile) || !file.exists(outfile)) {
		
		## this is a new PAML run
		## change working directory, since PAML clutters it up
		wd.old <- getwd()
		wd <- tempdir()
		setwd(wd)
		
		## prepare input files
		infile <- paste0(tempfile(), ".in")
		if (!is.character(aln)) {
			if (!inherits(aln, "DNAbin"))
				aln <- as.DNAbin(aln)
			aln <- toupper(makeLabel(aln))
			.write.paml.input(aln, infile, colsep = " ", indent = 10)
		}
		else {
			infile <- normalizePath(aln)
		}
		if (!is.character(tree)) {
			treefile <- paste0(infile, ".tree")
			write.tree(tree, treefile)
		}
		else {
			treefile <- normalizePath(tree)
		}
		
		ctl <- paste0(infile, ".ctl")
		outfile <- paste0(infile, ".out")
		
		## make contents of control file
		# i/o basics
		ll <- c(paste("seqfile =", infile),
				paste("treefile =", treefile),
				paste("outfile =", outfile),
				paste("noisy = 9"),
				paste("verbose = 1"),
				paste("runmode = 0"))
		# model specification
		ll <- c(ll,
				paste("seqtype =", seqtype),
				paste("clock =", clock),
				paste("model =", model),
				paste("NSsites =", ns.sites),
				paste("omega =", omega),
				paste("fix_omega =", fix.omega))
		
		## write control file
		writeLines(ll, ctl)
		
		## run PAML
		message("Running PAML::codeml in directory <", wd,">")
		message("\tcontrol file:\t", basename(ctl))
		message("\tsequence file:\t", basename(infile))
		message("\ttree file:\t", basename(treefile))
		message("\toutput file:\t", basename(outfile))
		cmd <- paste("codeml", ctl)
		stuff <- system(cmd, intern = !verbose)
		
		## return to previous working directory
		setwd(wd.old)
		
	}
	else {
		## extract stuff from an existing run
		infile <- NA
		ctl <- NA
		treefile <- NA
	}
	
	## read result (2lnL, df, ...)
	output <- readLines(outfile)
	x <- strsplit(output[ grepl("lnL", output) ], "[\\)\\:]*\\s+", perl = TRUE)
	rez <- list()
	rez$alignment <- infile
	rez$tree <- treefile
	rez$control <- ctl
	rez$output <- outfile
	rez$df <- as.numeric(x[[1]][4])
	rez$lnL <- as.numeric(x[[1]][5])
	
	## done
	return(rez)
	
}