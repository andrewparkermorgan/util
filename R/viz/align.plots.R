## Author:	Baptiste Augiue (baptiste.auguie@gmail.com)
## Date:		22 March 2013
## https://gist.github.com/baptiste/2877821

require(grid)
require(ggplot2)

set_panel_size <- function(p=NULL, g=ggplotGrob(p), width=unit(3, "cm"), height=unit(3, "cm")){
	panel_index_w<- g$layout$l[g$layout$name=="panel"]
	panel_index_h<- g$layout$t[g$layout$name=="panel"]
	g$widths[[panel_index_w]] <- width
	g$heights[[panel_index_h]] <- height
	class(g) <- c("fixed", class(g), "ggplot")
	g
}

print.fixed <- function(x) grid.draw(x)

left_width <- function(g){
	axis_l_index <- g$layout$r[g$layout$name=="axis-l"]
	ylab_index <- g$layout$r[g$layout$name=="ylab"]
	g$widths[[axis_l_index]] + g$widths[[ylab_index]]
}

full_width <- function(g){
	sum(g$widths)
}


align.plots <- function(..., width=unit(3, "cm"), height=unit(1, "null")){
	
	pl <- list(...)
	gl <- lapply(pl, set_panel_size, width=width, height=height)
	
	left <- lapply(gl, left_width)
	max_left <- max(do.call(unit.c, left))
	
	widths <- lapply(gl, full_width)
	max_width <- max(do.call(unit.c, widths))
	
	lay <- grid.layout(nrow=length(gl), ncol=1)
	vp <- viewport(layout=lay)
	pushViewport(vp)
	
	for(ii in seq_along(gl)){
		pushViewport(viewport(layout.pos.row=ii))
		pushViewport(viewport(x=unit(0.5, "npc") - 0.5*max_width + max_left - left[[ii]],
													just="left", width=widths[[ii]]))
		grid.draw(gl[[ii]])
		upViewport(2)
	}
	upViewport()

}