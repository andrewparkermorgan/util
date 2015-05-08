## some functions for doing common stuff with strings

## in the style of python's string.format()
formatstr <- function(x, ...) {
	
	args <- list(...)
	hasnames <- !is.null(names(args))
	
	for (i in seq_along(args)) {
		if (hasnames)
			x <- gsub(paste0("{", names(args)[i], "}"), args[i], x, fixed = TRUE)
		else
			x <- sub("{}", args[i], x, fixed = TRUE)
	}
	
	return(x)
	
}