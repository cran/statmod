mlreg.fit <-
#  Fit a linear model by maximum likelihood.
#  Distribution is one of "weibull", "exponential", "gaussian", "logistic",
#  "lognormal" or "loglogistic".
#  See survreg and survreg.object in the survival library for further documentation.
#	Gordon Smyth, smyth@wehi.edu.au
#  21 Nov 2001

function(X, y, weights=NULL, dist="logistic", init=NULL, scale=NULL) {
	require("survival")
	dist <- match.arg(dist,c("weibull","exponential","gaussian","logistic","lognormal","loglogistic"))

#	Remove rows with missing values
	if(any(is.na(X)) || any(is.na(y))) {
		anyna <- apply(is.na(X),1,any) | is.na(y)
		omit <- seq(along=y)[anyna]
		X <- X[-omit,]
		y <- y[-omit]
	}

	Y <- cbind(y,1)
	if(is.null(scale)) scale <- 0
	dist <- survreg.distributions[[dist]]
	controlvals <- survreg.control()
	offset <- rep(0,length(y))
	survival:::survreg.fit(X, Y, weights=weights, dist=dist, init=init, scale=scale, controlvals=controlvals, offset=offset, nstrat=1, strata=0, parms=NULL)
}

mlreg.fit.zero <- 
#  Fit a zero mean model by maximum likelihood (to get scale estimate and likelihood value).
#  Distribution is one of "weibull", "exponential", "gaussian", "logistic",
#  "lognormal" or "loglogistic".
#  See survreg and survreg.object in the survival library for further documentation.
#	Gordon Smyth, smyth@wehi.edu.au
#  21 Nov 2001.  Last modified 10 Feb 2004.

function(y, weights=NULL, dist="logistic", init=NULL, scale=NULL) {
	require("survival")
	dist <- match.arg(dist,c("weibull","exponential","gaussian","logistic","lognormal","loglogistic"))

#	Remove missing observations
	if(any(is.na(y))) y <- y[!is.na(y)]
	n <- length(y)

#	survreg can't handle zero means, so double-up y to ensure estimated mean is zero
	y <- c(y,-y)

	Y <- cbind(y,1)
	X <- matrix(1,length(y),1)
	if(is.null(scale)) scale <- 0
	dist <- survreg.distributions[[dist]]
	controlvals <- survreg.control()
	offset <- rep(0,length(y))
	out <- survival:::survreg.fit(X, Y, weights=weights, dist=dist, init=init, scale=scale, controlvals=controlvals, offset=offset, nstrat=1, strata=0, parms=NULL)

#	Now restore the output as for the original number of observations
	out$coefficients <- out$coefficients[2]
	out$icoef <- out$coefficients
	out$var <- 2 * out$var[2,2,drop=F]
	out$loglik <- out$loglik / 2
	out$linear.predictors <- out$linear.predictors[1:n]
	out$df <- 1
	out
}
