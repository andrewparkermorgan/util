MM.intensity.plot <- function(ids = NULL, sample.group = NULL, geno.class = NULL, target = NULL, markers = NULL, by = "name", ...) {
	
	require(GenomicRanges)
	require(ggplot2)
	require(plyr)
	
	## fetch samples from database
	samples <- fetch.samples(ids = ids, group = sample.group, by = by, ...)
	
	## fetch intensities in target region (will take a while if db is large)
	chr <- NULL
	start <- NULL
	end <- NULL
	if (inherits(target, "GRanges") & length(target)) {
		chr <- gsub("^chr","", seqnames(target)[1])
		start <- start(target)[1]
		end <- end(target)[1]
	}
	intens <- fetch.intensities(samples$id, chr = chr, start = start, end = end, by = "id", ...)
	if (!is.null(markers))
		intens <- subset(intens, marker %in% markers)
	
	## summarize intensities
	sums <- ddply(intens, .(sid), summarize, sum.intens = sum(x+y))
	sums <- merge(sums, samples[ ,c("id","name","batch") ], by.x = "sid", by.y = "id")
	sums$geno.class <- geno.class[ as.character(samples$name) ]
	
	## make plot (which is this function's return value)
	p <- ggplot(sums, aes(x = sum.intens))
	
	if (!is.null(geno.class))
		p <- p + geom_point(aes(y = reorder.by.group(sid, geno.class, sum.intens), colour = geno.class))
	else
		p <- p + geom_point(aes(y = reorder(sid, sum.intens)))

	p <- p + theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank()) +
		ylab("samples (one per row)\n") + xlab("\nsum-intensity") +
		ggtitle.n(paste("sum-intensity in", .region.to.str(target)), nrow(sums))
	
}

.region.to.str <- function(gr, ...) {
	
	stopifnot(inherits(gr, "GRanges") & length(gr))	
	paste(seqnames(gr)[1], ":", round(start(gr)[1]/1e6, 2), "-", round(end(gr)[1]/1e6, 2), "Mbp")

}

reorder.by.group <- function(x, groups, values, ...) {
	
	if (!is.factor(groups))
		groups <- factor(groups)
	
	rez <- lapply(levels(groups), function(g) {
		i <- (groups == g)
		reorder(x[i], values[i], ...)
	})
	
	lev <- do.call("c", lapply(rez, levels))
	# print(lev)
	factor(x, levels = lev)
	
}