library(Matrix)
library(plyr)
library(ggplot2)
library(RColorBrewer)

source("~/lib/util/R/intervals.R")
source("~/lib/util/R/mouse.R")
source("~/lib/util/R/viz/genome.r")

load("~/db/subspecies/sso.Rdata")
load("~/db/populations/CC/ss")
strains <- read.table("~/tmp/strains.txt", header = TRUE, sep = "\t")

## retrieve row-col indices for nonzero elements of a dsCMatrix
get.sparse.indices <- function(x, ...) {
	
	if (!inherits(x, "dsCMatrix"))
		stop("I was expecting an object of class 'dsCMatrix.'")
	
	first <- diff(x@p)
	rr <- x@i
	cc <- do.call("c", sapply(1:ncol(x), function(i) rep(i, first[i])) )
	return( cbind(rr,cc) )
	
}

bins <- make.windows(mm9, 0.5e6)
bins <- bins[ start(bins) > 3e6 & !(seqnames(bins) %in% c("chrY","chrM")) ]
sso.binned <- llply(sso.reduced, function(x) {
	tryCatch({
		project.onto(x, bins)
	}, error = function(e) {
		GRanges()
	})
}, .progress = "text") # takes forever
pairs <- combn(c("dom","mus","cas"), 2)

rez <- alply(pairs, 2, function(p) {

	#as.character(subset(strains, group == "classical")$strain)
	#stacked <- llply(as.character(subset(strains, group == "classical")$strain), function(s) {
	stacked <- llply(names(sso.cc), function(s) {
	
		## project subspecies blocks onto bins
		#z <- sso.binned[[s]]
		z <- sso.cc[[s]][ as.vector(seqnames(sso.cc[[s]])) == "chr1" ]
		
		## compute incidence matrix
		i <- as.integer(values(z)$ss == p[1])
		i[ is.na(i) ] <- 0
		j <- as.integer(values(z)$ss == p[2])
		j[ is.na(j) ] <- 0
		is <- as(i, "sparseMatrix")
		js <-  as(j, "sparseMatrix")
		X <- tcrossprod(is,js)
		Xp <- tcrossprod(js,is)
		return( forceSymmetric(X+Xp) )
		
	}, .progress = "text")
	
	return( Reduce("+", stacked) )
	
})

allpairs <- data.frame()
pairlabels <- apply(t(pairs), 1, paste, collapse = "/")
for (i in seq_along(rez)) {
	
	df <- as.data.frame(cbind(get.sparse.indices(rez[[i]]), rez[[i]]@x))
	colnames(df) <- c("x","y","z")
	df$chr <- as.vector(seqnames(bins))[ df$x+1 ]
	df$start <- start(bins)[ df$x+1 ]
	df$end <- start(bins)[ df$x+1 ]
	allpairs <- rbind(allpairs, transform(df, pair = pairlabels[i]))
	
}

chroms <- ddply(allpairs, .(chr), summarize, first = min(x), last = max(x), mid = floor(mean(range(x))))
chroms$chr <- factor(chroms$chr, levels = seqnames(mm9)[1:20])
chroms <- chroms[ order(chroms$chr), ]

## plot one chromosome at a time (screen okay)
ggplot(subset(allpairs, chr == "chr1")) +
	geom_tile(aes(x = x, y = y, fill = pair, alpha = z)) +
	geom_abline(xintercept = 0, yintercept = 0, slope = 1, lty = "dashed", col = "grey") +
	#geom_rect(data = chroms, aes(ymin = -100, ymax = 0, xmin = first-1, xmax = last, alpha = chr)) +
	#scale_alpha_manual(values = rep_len(c(0.2,1), nlevels(chroms$chr))) +
	scale_alpha_continuous(range = c(0.3,1)) +
	scale_fill_manual("subspecies pair", values = brewer.pal(3, "Set1")[c(1,3,2)]) +
	scale_x_continuous(limits = c(0, max(allpairs$x)), labels = function(x) paste(round(x*0.5, 1))) +
	scale_y_continuous(limits = c(0, max(allpairs$x)), labels = function(x) paste(round(x*0.5, 1))) +
	facet_grid(pair ~ .) +
	coord_equal() +
	#scale_x_continuous(breaks = chroms$mid, labels = gsub("^chr","", chroms$chr)) +
	#scale_y_continuous(limits = c(-100, max(chroms$last)+100)) +
	guides(alpha = FALSE) +
	#ggtitle("pairwise subspecies combinations in 100 classical inbred strains\n") +
	ggtitle("pairwise subspecies combinations in 69 CC lines\n") +
	xlab("\nposition (Mbp)") + ylab("position (Mbp)\n") +
	theme_bw()

## plot genome-wide (need to send straight to file, to avoid overwhelming screen device)
png("~/tmp/heatmap.png", width = 1200, height = 1200)
ymax <- max(allpairs$y) + 0.05*max(allpairs$y)
ymin <- -0.05*max(allpairs$y)
ggplot(allpairs) +
	geom_tile(aes(x = x, y = y, fill = pair, alpha = z)) +
	geom_segment(data = chroms, aes(y = ymax, yend = ymax, x = first+1, xend = last, colour = chr)) +
	geom_segment(data = chroms, aes(x = ymin, xend = ymin, y = first+1, yend = last, colour = chr)) +
	scale_alpha_continuous(range = c(0.3,1)) +
	scale_colour_manual(values = rep_len(c("grey80","grey30"), nlevels(chroms$chr))) +
	scale_fill_manual("subspecies pair", values = brewer.pal(3, "Set1")[c(1,3,2)]) +
	scale_x_continuous(limits = c(ymin, ymax), breaks = chroms$mid, labels = gsub("^chr","", chroms$chr)) +
	scale_y_continuous(limits = c(ymin, ymax), breaks = chroms$mid, labels = gsub("^chr","", chroms$chr)) +
	geom_abline(xintercept = 0, yintercept = 0, slope = 1, lty = "dashed", col = "grey") +
	coord_equal() +
	guides(alpha = FALSE, colour = FALSE) +
	ggtitle("pairwise subspecies combinations in 100 classical inbred strains\n") +
	theme_axesonly() + theme(axis.title = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
dev.off()