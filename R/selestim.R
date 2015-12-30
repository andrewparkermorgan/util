## selestim.R
# helper functions for using the SelEstim software (http://www1.montpellier.inra.fr/CBGP/software/selestim/)

make.selestim <- function(gty, path = tempfile(pattern = "selestim", fileext = ".dat"), ...) {
	
	if (!inherits(gty, "genotypes"))
		stop("Please supply an object of class 'genotypes'.")
	
	## get allele counts (major+minor) by marker, within populations
	message("Calculating allele counts...")
	f <- genoapply(gty, 2, .(fid), function(x) {
		mac <- freq(x, counts = TRUE)
		n <- 2*ncol(x)
		cbind(mac, n-mac)
		})
	ff <- do.call("cbind", f)
	#print(head(ff))
	
	## write output
	message("Writing to file <", path, ">...")
	outfile <- file(path, open = "w")
	npop <- ncol(ff)/2
	nsnps <- nrow(ff)
	cat(as.character(npop), "\n", file = outfile)
	cat(as.character(nsnps), "\n", file = outfile)
	write.table(ff, outfile, append = TRUE, sep = " ",
				quote = FALSE, row.names = FALSE, col.names = FALSE)
	close(outfile)
	
	message("Done.")
	invisible(path)
	
}