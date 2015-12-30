## layout.R
## functions for arranging ggplot objects on a page

library(grid)
library(gridExtra)

grid.arrange.list <- function(pl, ...) {
	
	call <- sapply(1:length(pl), function(i) paste0("pl[[",i,"]]"))
	call <- paste(call, collapse = ",")
	call <- paste0("grid.arrange(", call, ", ...)")
	eval(parse(text = call))

}
