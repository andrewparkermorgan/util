library(proto)
library(grid)
library(ggplot2)



geom_track <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity", space = 0.1, ...) { 
	GeomTrack$new(mapping = mapping, data = data, stat = stat, position = position, space = space, ...)
}

GeomTrack <- proto(ggplot2:::Geom, {
	
	draw <- function(., data, scales, coordinates, space = 0.1, ...) 
	{
		
		data$xmin <- data$start
		data$xmax <- data$end
		data$ymin <- 1-space
		data$ymax <- 0
		if ("group" %in% names(data)) {
			data$offset <- as.numeric(factor(data$group))
		}
		data <- transform(data, ymin = offset-1, ymax = offset-space)
		ggplot2:::GeomRect$draw_groups(data, scales, coordinates)
	
	}
	
	# draw_legend <- function(., data, ...) {
		# show a grob in the key that represents the data
	# }
	
	objname <- "track" # name of the geom in lowercase. Must correspond to GeomField.
	desc <- "Multi-track underlay of genomic intervals"
	
	default_stat <- function(.) ggplot2:::StatIdentity
	required_aes <- c("start","end")
	default_aes <- function(.) ggplot2:::aes(colour = "black") #ggplot2:::GeomRect$default_aes()
	guide_geom <- function(.) "rect"
	
})