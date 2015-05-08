## --- hmm.R --- ##
##	Date: 10 March 2015
##	Purpose: Finally implement a general HMM for myself in pure R, allowing time-varying transition and emission probabilities.

rescale <- function(x, y = 1) {
	if (sum(x) > 0)
		y*x/sum(x)
	else
		rep(0, length(x))
}

## initialize an HMM object with known transition and emission probabilities
init.hmm <- function(x, tprob, emit, ...) {
	
	## set length of observed sequence
	if (!length(x))
		stop("Observed sequence must have strictly positive length.")
	else if (length(x) == 1)
		nobs <- as.integer(ceiling(abs(x)))
	else
		nobs <- length(x)
	
	## check validity of transition matrix
	if (length(dim(tprob)) < 2 | !(is.matrix(tprob)|is.array(tprob)) | !is.numeric(tprob))
		stop("Transition probabilities are likely misspecified; please supply a numeric matrix with >=2 dimensions.")
	else if (length(dim(tprob)) < 3) {
		## if transition probs are not time-varying, pretend like they are
		.tprob <- array(dim = c(nobs, dim(tprob)))
		for (i in 1:nobs)
			.tprob[i,,] <- tprob
		tprob <- .tprob
	}
	
	## normalize transition matrix to columns and rows sum to unity
	for (i in 1:nobs) {
		if (!identical(tprob[i,,], t(tprob[i,,])))
			stop("Oops -- transition matrix should be symmetric.")
		tprob[i,,] <- tprob[i,,]/rowSums(tprob[i,,])
		if ( any(colSums(tprob[i,,]) != 1) | any(rowSums(tprob[i,,]) != 1) )
			stop("Normalization of transition matrix failed.")
	}
	
	if (length(dim(emit)) < 2 | !is.matrix(emit) | !is.numeric(emit))
		stop("Emission probabilities are likely misspecified; please supply a numeric matrix with >= 2 dimensions.")
	else if (length(dim(emit)) < 3) {
		## if transition probs are not time-varying, pretend like they are
		.emit <- array(dim = c(nobs, dim(emit)))
		for (i in 1:nobs)
			.emit[i,,] <- emit
		emit <- .emit
	}
	
	## normalize emission probabilities so that rows sum to unity
	for (i in 1:nobs) {
		emit[i,,] <- emit[i,,]/rowSums(emit[i,,])
		if ( any(rowSums(tprob[i,,]) != 1) )
			stop("Normalization of emission matrix failed.")
	}
	
	#print(any(is.na(tprob)))
	#print(any(is.na(emit)))
	
	nstates <- dim(emit)[2]
	nsymb <- dim(emit)[3]
	obs <- NULL
	if (length(x) > 1)
		obs <- x
	
	if (length(x) > 1)
		if (length(unique(x)) != nsymb)
			warning(paste0("Number of observed symbols (",length(unique(x)),
						   ") is different than corresponding dimension of emission matrix (",nsymb,"). ",
						   "This could be okay but you should check the input."))
	
	hmm <- list(nstates = nstates, nsymb = nsymb, nobs = nobs,
				tprob = tprob, emit = emit, obs = obs)
	class(hmm) <- c("hmm", class(hmm))
	return(hmm)
	
}

## get posterior decoding, given transition and emission probs in <hmm>
viterbi <- function(hmm, x = NULL, progress = TRUE, ...) {

	if (!inherits(hmm, "hmm"))
		warning("Are you sure this is an object of class 'hmm'?  Proceeding with skepticism...")
	
	## initalize observed state sequence
	if (is.null(x))
		x <- hmm$obs
	else
		if (length(x) != hmm$nobs)
			stop("Observation sequence provided doesn't match the one with which the model was intialized.")
	x <- as.numeric(factor(x, nmax = hmm$nsymb))
	
	log1p <- function(x) log(1+x)
	
	## from Christinanini & Hahn (2007) pp75-77
	V <- matrix(0, nrow = hmm$nstates, ncol = hmm$nobs+1)
	tprob <- hmm$tprob
	emit <- hmm$emit
	V[,1] <- log1p(1/nrow(V))
	tb <- matrix(0, nrow = hmm$nstates, ncol = hmm$nobs)
	if (progress)
		pb <- txtProgressBar(min = 0, max = ncol(V)-1, style = 3)
	for (i in 2:ncol(V)) {
		for (k in 1:nrow(V)) {
			mu <- V[ ,i-1 ] + log1p(tprob[ i-1,,k ])
			#best <- which.max(mu)
			#tb[ k,i-1 ] <- best
			V[ k,i ] <- max(mu) + log1p(emit[ i-1,k,x[i-1] ])
			last <- which.max(V[,i])
		}
		if (progress)
			setTxtProgressBar(pb, i-1)
	}
	
	#print(tb)
	#print(last)
	#decoded <- rep(0, hmm$nobs)
	#decoded[ hmm$nobs ] <- last
	#for (i in ncol(V):2) {
	#	decoded[i] <- tb[ last, i-1 ]
	#	last <- tb[ last,i-1 ]
	#}
	decoded <- apply(V, 2, which.max)
	
	hmm$decoded <- decoded[-1]
	hmm$encoded <- x
	
	return(hmm)
	
}

simulate.hmm <- function(hmm, n, init = NULL, ...) {
	
	if (!inherits(hmm, "hmm"))
		warning("Are you sure this is an object of class 'hmm'?  Proceeding with skepticism...")
	
	if (n < 1)
		stop("Can only simulate sequences of length >= 1.")
	if (n > hmm$nobs)
		stop("Can only simulate sequences of length <= nobs.  Initialize a longer HMM if you need more.")
	
	tprob <- hmm$tprob
	emit <- hmm$emit
	
	## generate hidden state sequence
	if (is.null(init))
		init <- rep(1/hmm$nstates, hmm$nstates)
	states <- rep(0, n)
	states[1] <- apply(rmultinom(1, 1, init), 2, which.max)
	for (j in 2:n) {
		states[j] <- apply(rmultinom(1, 1, tprob[j-1,states[j-1],]), 2, which.max)
	}
	
	## simulate observed states conditional on hidden states
	obs <- rep(0, n)
	for (i in 1:n) {
		obs[i] <- apply(rmultinom(1, 1, emit[ i,states[i], ]), 2, which.max)
	}
	
	return( list(hidden = states, obs = obs) )
	
}