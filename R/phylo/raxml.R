## some functions to control RAxML

raxml.dnaml <- function(aln, model = "GTRGAMMA", ...) {
	
	rez <- .run.raxml(aln, model = model, ...)
	return(rez)
	
}

raxml.multistate <- function(aln, model = "MULTIGAMMA", ...) {
	
	rez <- .run.raxml(aln, model = model, ...)
	return(rez)
	
}

raxml.parsimony <- function(aln, ...) {
	
	rez <- .run.raxml(aln, model = NULL, outprefix = "RAxML_parsimonyTree", ...)
	return(rez)
	
}

.write.fasta <- function(aln, outfile, clean.names = FALSE, ...) {
	
	if (clean.names)
		aln$nam <- make.names(aln$nam)
	
	ff <- file(outfile, open = "w")
	for (i in seq_len(aln$nb)) {
		writeLines(paste0(">", aln$nam[i]), ff)
		writeLines(paste0(aln$seq[i]), ff)
	}
	close(ff)
	
}

.run.raxml <- function(aln, outprefix = "RAxML_bipartitions", path = "raxml", model = NULL, nreps = 100,
					   boot.seed = 12345, random.seed = 12345, outgroup = NULL,
					   flags = "-T 2", ...) {
	
	ff <- tempfile()
	wd <- dirname(ff)
	infile <- paste0(ff, ".in")
	outfile <- paste0(ff, ".out")
	
	
	cmd <- paste(path, "-n", basename(outfile), "-s", (infile),
				 "-w", wd, "-p", random.seed)
	if (!any(is.null(model)))
		cmd <- paste(cmd, "-f a", "-m", model,
					 "-x", boot.seed, "-#", nreps)
	else
		cmd <- paste(cmd, "-m GTRGAMMA", "-y")
	
	old.names <- aln$nam
	new.names <- make.names(aln$nam)
	renamer <- setNames( old.names, new.names )
	aln$nam <- new.names
	
	.write.fasta(aln, infile, clean.names = TRUE)
	if (!is.null(outgroup))
		cmd <- paste(cmd, "-o", paste(outgroup, collapse = ","))
	
	cmd <- paste(cmd, flags)
	print(cmd)
	rez <- system(cmd, intern = FALSE)
	
	if (!rez) {
		treefile <- file.path(wd, paste0(outprefix, ".", basename(outfile)))
		tree <- read.tree(treefile)
		tree$tip.label <- renamer[ tree$tip.label ]
		attr(tree, "alignment") <- infile
		attr(tree, "file") <- treefile
		return(tree)
	}
	else {
		return(FALSE)
	}
	
}