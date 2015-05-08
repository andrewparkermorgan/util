## --- treemix.R ---
## 1 March 2015
## utilities for handling i/o for TreeMix runs
## NB: requires functions in <plink.R>

## compute allele frequencies with plink, convert to treemix format
make.treemix <- function(prefix, outfile = paste0(as.character(prefix), ".treemix.gz"), freq = NULL, ...) {
	
	require(reshape2)
	
	if (1) {
		cmd <- paste0("--freq counts --family --out ", prefix)
		freqfile <- paste0(prefix, ".frq.strat")
		success <- .plink.command(prefix, cmd, list(freqfile))
		if (!success)
			stop("Attempt to compute allele frequencies failed.")	
	}
	
	if (is.null(freq))
		freq <- read.table(freqfile, header = TRUE)
	totals <- acast(freq, SNP ~ CLST, value.var = "NCHROBS")
	minors <- acast(freq, SNP ~ CLST, value.var = "MAC")
	majors <- totals - minors
	
	out <- gzfile(outfile)
	open(out, "w")
	on.exit(close(out))
	
	writeLines(paste(colnames(totals), collapse = " "), out)
	nr <- nrow(totals)
	pb <- txtProgressBar(min = 0, max = nr, style = 3)
	for (i in 1:nr) {
		line <- paste(majors[i,], minors[i, ], sep = ",", collapse = " ")
		writeLines(line, out)
		setTxtProgressBar(pb, i)
	}
	
}

## read results from a TreeMix run and compute pseudo-R^2
read.treemix <- function(stem, ...) {
	
	require(ape)
	
	ff <- sapply(c("cov","modelcov","treeout"), function(f) paste0(stem, ".", f, ".gz"))
	if (!all(sapply(ff, file.exists))) {
		stop(paste0("Can't find all outputs from TreeMix run with prefix '",stem,"'"))
	}
	
	cov <- as.matrix( read.table(gzfile(ff[1]), as.is = TRUE, head = TRUE, quote = "", comment.char = "") )
	mod <- as.matrix( read.table(gzfile(ff[2]), as.is = TRUE, head = TRUE, quote = "", comment.char = "") )
	resid <- mod - cov
	i <- upper.tri(resid, diag = FALSE)
	sse <- sum( (resid[i] - mean(resid[i]))^2 )
	ssm <- sum( (mod[i] - mean(mod[i]))^2 )
	r2 <- 1 - sse/ssm
	
	tree <- read.tree(text = readLines(gzfile(ff[3]), n = 1))
	
	return( list(cov = cov, cov.est = mod, resid = resid,
				 sse = sse, ssm = ssm, r2 = r2,
				 tree = tree) )
	
}

## draw the characteristic TreeMix tree with migration edges
plot.treemix <- function(stem, o = NA, flip = vector(), draw = TRUE, branch.colour = "grey", 
						 plot.nodes = TRUE, plot.migration = TRUE, ...) {
	
	d <- paste0(stem, ".vertices.gz")
	e <- paste0(stem, ".edges.gz")
	#se <- paste0(stem, ".covse.gz")
	d <- read.table(gzfile(d), as.is = TRUE, comment.char = "", quote = "")
	e <- read.table(gzfile(e), as.is  = TRUE, comment.char = "", quote = "")
	
	e[,3] <- e[,3]*e[,4]
	e[,3] <- e[,3]*e[,4]
	
	#se <- read.table(gzfile(se), as.is = TRUE, comment.char = "", quote = "")
	#m1 <- apply(se, 1, mean)
	#m <- mean(m1)
	for(i in seq_along(flip)) {
		d <- .flip.node(d, flip[i])
	}
	d$x <- "NA"
	d$y <- "NA"
	d$ymin <- "NA"
	d$ymax <- "NA"
	d$x <- as.numeric(d$x)
	d$y <- as.numeric(d$y)
	d$ymin <- as.numeric(d$ymin)
	d$ymax <- as.numeric(d$ymax)
	
	d <- .set.y.coords(d)
	d <- .set.x.coords(d, e)
	d <- .set.mig.coords(d, e)
	
	print(e)
	stuff <- .layout.treemix(d, e, o = o)
	
	tdepth <- max(stuff$tips$x, na.rm = TRUE)
	d$xo <- d$x + 0.01*tdepth
	stuff$tips$xo <- stuff$tips$x + 0.01*tdepth
	xlim <- c(0, max(stuff$tips$xo, na.rm = TRUE)*1.05)
	p <- ggplot(stuff$tips) +
		geom_segment(data = subset(stuff$edges, type == "NOT_MIG"),
					 aes(x = from.x, y = from.y, xend = to.x, yend = to.y),
					 colour = branch.colour) +
		geom_text(aes(x = xo, y = y, label = pop), hjust = 0) +
		scale_x_continuous(limits = xlim) +
		xlab("\ndrift parameter")

	if (plot.nodes)
		p <- p + geom_point(data = d, aes(x = x, y = y))
	
	if (plot.migration) {
		p <- p + geom_segment(data = subset(stuff$edges, type == "MIG"),
							  aes(x = from.x, y = from.y, xend = to.x, yend = to.y, colour = weight),
							  #arrow = arrow(length = unit(6, "points")),
							  alpha = 0.5) +
			scale_colour_gradient2("migration weight", high = "red", mid = "yellow")
	}
		
	
	if (draw)
		print(p)
	
	return( list(vertices = d, edges = e, plot = p) )
	
}


### modified from TreeMix source

.set.y.coords <- function(d) {
	
	i <- which(d[,3]=="ROOT")
	y <- d[i,8]/ (d[i,8]+d[i,10])
	d[i,]$y <- 1-y
	d[i,]$ymin <- 0
	d[i,]$ymax <- 1
	c1 <- d[i,7]
	c2 <- d[i,9]
	ni <- which(d[,1]==c1)
	ny <- d[ni,8]/ (d[ni,8]+d[ni,10])
	d[ni,]$ymin <- 1-y
	d[ni,]$ymax <- 1
	d[ni,]$y <- 1- ny*(y)
	
	ni <- which(d[,1]==c2)
	ny <- d[ni,8]/ (d[ni,8]+d[ni,10])
	d[ni,]$ymin <- 0
	d[ni,]$ymax <- 1-y
	d[ni,]$y <- (1-y)-ny*(1-y)
	
	for (j in 1:nrow(d)){
		d <- .set.y.coord(d, j)
	}	
	return(d)
}

.set.y.coord <- function(d, i) {
	
	index <- d[i,1]
	parent <- d[i,6]
	if (!is.na(d[i,]$y)){
		return(d)
	}
	tmp <- d[d[,1] == parent,]
	if ( is.na(tmp[1,]$y)){
		d <- .set.y.coord(d, which(d[,1]==parent))
		tmp <- d[d[,1]== parent,]
	}
	py <- tmp[1,]$y
	pymin <- tmp[1,]$ymin
	pymax <- tmp[1,]$ymax
	f <- d[i,8]/( d[i,8]+d[i,10])
	#print (paste(i, index, py, pymin, pymax, f))
	if (tmp[1,7] == index){
		d[i,]$ymin <- py
		d[i,]$ymax <- pymax
		d[i,]$y <- pymax-f*(pymax-py)
		if (d[i,5]== "TIP"){
			d[i,]$y <- (py+pymax)/2
		}
	}
	else{
		d[i,]$ymin <- pymin
		d[i,]$ymax <- py
		d[i,]$y <- py-f*(py-pymin)
		if (d[i,5]== "TIP"){
			d[i,]$y <- (pymin+py)/2
		}	
		
	}
	return(d)
}


.set.x.coords <- function(d, e) {
	
	i <- which(d[,3]=="ROOT")
	index <- d[i,1]
	d[i,]$x <- 0
	c1 <- d[i,7]
	c2 <- d[i,9]
	ni <- which(d[,1]==c1)
	tmpx <-  e[e[,1]==index & e[,2] == c1,3]
	if (length(tmpx) == 0){
		tmp <- e[e[,1] == index,]
		tmpc1 <- tmp[1,2]
		if ( d[d[,1]==tmpc1,4] != "MIG"){
			tmpc1 <- tmp[2,2]
		}
		tmpx <- .get.dist.to.nmig(d, e, index, tmpc1)
	}
	if(tmpx < 0){
		tmpx <- 0
	}
	d[ni,]$x <- tmpx
	
	ni <- which(d[,1]==c2)
	tmpx <-  e[e[,1]==index & e[,2] == c2,3]
	if (length(tmpx) == 0){
		tmp <- e[e[,1] == index,]
		tmpc2 <- tmp[2,2]
		if ( d[d[,1]==tmpc2,4] != "MIG"){
			tmpc2 <- tmp[1,2]
		}
		tmpx = .get.dist.to.nmig(d, e, index, tmpc2)
	}
	if(tmpx < 0){
		tmpx <- 0
	}
	d[ni,]$x <- tmpx
	
	for (j in 1:nrow(d)){
		d <- .set.x.coord(d, e, j)
	}
	return(d)

}


.set.x.coord <- function(d, e, i) {
	
	index <- d[i,1]
	parent <- d[i,6]
	if (!is.na(d[i,]$x)){
		return(d)
	}
	tmp <- d[d[,1] == parent,]
	if ( is.na(tmp[1,]$x)){
		d <- .set.x.coord(d, e, which(d[,1]==parent))
		tmp <- d[d[,1]== parent,]
	}
	#print (paste(parent, index))
	tmpx <- e[e[,1]==parent & e[,2] == index,3]
	if (length(tmpx) == 0){
		tmp2 <- e[e[,1] == parent,]
		tmpc2 <- tmp2[2,2]
		#print
		if ( d[d[,1]==tmpc2,4] != "MIG"){
			tmpc2 <- tmp2[1,2]
		}
		tmpx <- .get.dist.to.nmig(d, e, parent, tmpc2)
	}
	if(tmpx < 0){
		tmpx <- 0
	}
	d[i,]$x <- tmp[1,]$x+ tmpx
	return(d)
	
}

.layout.treemix <- function(d, e, o = NA, ...) {
	
	## containers for return values
	ne <- nrow(e)
	from.x <- numeric(ne)
	from.y <- numeric(ne)
	to.x <- numeric(ne)
	to.y <- numeric(ne)
	etype <- character(ne)
	weight <- numeric(ne)
	
	## loop over edges
	for(i in 1:nrow(e)) {
		v1 <- d[d[,1] == e[i,1],]
		v2 <- d[d[,1] == e[i,2],]
		from.x[i] <- v1[1,]$x
		from.y[i] <- v1[1,]$y
		to.x[i] <- v2[1,]$x
		to.y[i] <- v2[1,]$y
		etype[i] <- e[i,5]
		weight[i] <- e[i,4]
	}
	
	edges <- data.frame(from.x = from.x, from.y = from.y, to.x = to.x, to.y = to.y,
						type = factor(etype), weight = weight)
	
	tmp <- d[d[,5] == "TIP",]
	tips <- with(tmp, data.frame(x = x, y = y, pop = V2))
	
	return( list(edges = edges, tips = tips) )
	
}

.set.mig.coords = function(d, e) {
	
	for (j in 1:nrow(d)){
		if (d[j,4] == "MIG"){
			p <- d[d[,1] == d[j,6],]
			c <- d[d[,1] == d[j,7],]
			tmpe <- e[e[,1] == d[j,1],]
			y1 <- p[1,]$y
			y2 <- c[1,]$y
			x1 <- p[1,]$x
			x2 <- c[1,]$x
			
			mf <- tmpe[1,6]	
			if (is.nan(mf)){
				mf <- 0
			}
			#d[j,]$y = (y1+y2)* mf
			#d[j,]$x = (x1+x2) *mf
			d[j,]$y <- y1+(y2-y1)* mf
			#print(paste(mf, x1, x2))
			d[j,]$x <- x1+(x2-x1) *mf
		}	
		
	}
	
	return(d)
	
}

.get.f = function(stem){
	d <- paste(stem, ".cov.gz", sep = "")
	d2 <- paste(stem, ".modelcov.gz", sep = "")
	d <- read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
	d2 <- read.table(gzfile(d2), as.is = T, comment.char = "", quote = "")
	d <- d[order(names(d)), order(names(d))]
	d2 <- d2[order(names(d2)), order(names(d2))]
	tmpcf <- vector()
	tmpmcf <- vector()
	for (j in 1:nrow(d)){
		for (k in (j+1):nrow(d)){
			tmpcf <- append(tmpcf, d[j,k])
			tmpmcf <- append(tmpmcf, d[j,k] - d2[j,k])
		}
	}
	tmpv <- var(tmpmcf)/var(tmpcf)
	return(1-tmpv)
	
}

.get.dist.to.nmig <- function(d, e, n1, n2) {
	toreturn <- e[e[,1] == n1 & e[,2] == n2,3]
	#print(toreturn)
	while ( d[d[,1] == n2,4] == "MIG"){
		tmp <- e[e[,1] == n2 & e[,5] == "NOT_MIG",]
		toreturn <- toreturn+tmp[1,3]
		n2 <- tmp[1,2]
	}
	return(toreturn)
}

.flip.node = function(d, n){
	
	i <- which(d[,1] == n)
	t1 <- d[i,7]
	t2 <- d[i,8]
	d[i,7] <- d[i,9]
	d[i,8] <- d[i,10]
	d[i,9] <- t1
	d[i,10] <- t2
	return(d)
	
}