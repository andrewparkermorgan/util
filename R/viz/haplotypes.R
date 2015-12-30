# haplotypes.R
# plotting functions for haplotype reconstructions for CC/DO mice
# assumes argyle loaded, and mouse.R sourced

plot.hapfile <- function(df, ...) {
	
	df$phase <- factor(df$phase)
	df$updown <- ifelse(df$phase == levels(df$phase)[1], -1, 1)
	df$chromosome <- factor(df$chromosome, levels = seqnames(mm9))
	ggplot(df) +
		geom_rect(aes(xmin = start, xmax = end, ymin = updown, ymax = 0, fill = strain)) +
		facet_grid(chromosome ~ ., labeller = function(k,v) gsub("^chr","", v)) +
		scale_fill_manual(values = CC.COLORS) +
		scale_x_genome() +
		theme_gbrowse()
	
}