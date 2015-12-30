## admixture.R
# functions for handling output from ADMIXTURE (Alexander DH 2009 Genome Res)

library(reshape2)

#' Read ancestry proportions (`Q matrix`) from ADMIXTURE run
#' 
#' @details Assumes that estimation was performed from a PLINK *.bed file,
#' 	and that the corresponding *.fam file resides in same directory as output.
#' 
read.Q.matrix <- function(infile, ...) {
	
	if (grepl("\\.\\d+\\.Q$", infile)) {
		fam.path <- gsub("\\.\\d+\\.Q$",".fam", infile)
	}
	else {
		fam.path <- paste0(infile, ".fam")
		infile <- paste0(infile, ".Q")
	}
	
	message("Reading ancestry proportions from <", infile, ">...")
	Q <- read.table(infile, header = FALSE)
	
	message("Reading sample metadata from <", fam.path, ">...")
	fam <- argyle::read.fam(fam.path)
	
	rownames(Q) <- rownames(fam)
	colnames(Q) <- paste0("pop", seq_len(ncol(Q)))
	Q <- cbind(iid = rownames(Q), Q)
	
	class(Q) <- c("admixture", class(Q))
	return(Q)
	
}

#' Make standard stacked-bars plot of ancestry proportions.
#' 
#' @param cluster attempt to sort samples along axis by similarity
#' @param label show individual IDs along axis
#' 
plot.admixture <- function(Q, cluster = FALSE, label = FALSE, meta = NULL, ...) {
	
	if (cluster) {
		message("Attempting to sort samples attractively...")
		#d <- stats::dist(as.matrix(Q[,-1]), "minkowski")
		#cl <- hclust(d)
		#Q$iid <- reorder(factor(Q$iid), cl$order)
		Q <- Q[ do.call(order, as.list(Q[,-1])), ]
	}
	
	qq <- melt(Q, id.var = "iid")
	if (!is.null(meta)) {
		qq <- merge(qq, meta, by.x = "iid", by.y = "iid", all.x = TRUE)
	}
	
	qq <- plyr::arrange(qq, variable)
	p0 <- ggplot(qq) +
		geom_bar(aes(x = iid, y = value, fill = variable), stat = "identity") +
		theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
	
	if (!label)
		p0 <- p0 + theme(axis.ticks.x = element_blank(),
						 axis.text.x = element_blank())
	
	return(p0)
	
}