## genotype-by-taqman: assign delta Ct values in <x> to groups (by default, allow 2 clases)
## (really just a wrapper around a simple 1d Gaussian-mixture classifier)
classify.taqman <- function(x, clusters = 2, ...) {
	
	require(mclust)
	
	clusters <- mclustBIC(x, G = clusters, ...)
	rez <- summary(clusters, x)
	
	return( data.frame(class = rez$classification, prob = 1-rez$uncertainty) )
	
}

predict.copies <- function(normed, controls, ...) {
	
	require(lme4)
	
	m0 <- lm(log(copies) ~ resid, controls)
	
	b <- coef(m0)[1]
	m <- coef(m0)[2]
	ci <- confint(m0)
	
	cn <- with(normed, b + m*resid)
	lci <- with(normed, ci[ 1,1 ] + ci[ 2,1 ]*resid)
	uci <- with(normed, ci[ 1,2 ] + ci[ 2,2 ]*resid)
	
	return( data.frame(copies = exp(cn), copies.lo = exp(lci), copies.hi = exp(uci)) )
	
}