randomizedBlock <- function(formula, random, weights=NULL, fixed.estimates=TRUE, data=list(), subset=NULL, contrasts=NULL) {
#	REML for mixed linear models
#	Gordon Smyth, Walter and Eliza Hall Institute
#	28 Jan 2003

#	Extract model from formula
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	mf$fixed.estimates <- NULL
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	xvars <- as.character(attr(mt, "variables"))[-1]
	if((yvar <- attr(mt,"response")) > 0) xvars <- xvars[-yvar]
	xlev <- if(length(xvars) > 0) {
		xlev <- lapply(mf[xvars], levels)
		xlev[!sapply(xlev, is.null)]
	}
	y <- model.response(mf, "numeric")
	w <- model.weights(mf)
	x <- model.matrix(mt, mf, contrasts)
	random <- mf[["(random)"]]

#	Missing values not allowed
	if(any(is.na(y)) || any(is.na(x)) || any(is.na(random))) stop("Missing values not allowed")
	if(!is.null(weights)) if(any(is.na(weights))) stop("Missing values not allowed")

#	Design matrix for random effects
	lev <- unique.default(random)
	z <- 0 + (matrix(random,length(random),length(lev)) == t(matrix(lev,length(lev),length(random))))

	out <- randomizedBlockFit(y,x,z,w=w,fixed.estimates=fixed.estimates)
	if(fixed.estimates==TRUE) class(out) <- "lm"
	out
}

randomizedBlockFit <- function(y,X,Z,w=NULL,fixed.estimates=TRUE) {
#	Restricted maximum likelihood estimation for mixed linear models.
#	Fits the model  Y = X*BETA + Z*U + E  where BETA is fixed
#	and U is random.
#
#	GAMMA holds the variance components.  The errors E and
#	random effects U are assumed to have covariance matrices
#	EYE*GAMMA(1) and EYE*GAMMA(2) respectively.

#	Gordon Smyth, Walter and Eliza Hall Institute
#	Matlab version 19 Feb 94.  Converted to R, 28 Jan 2003.

#  Prior weights
if(!is.null(w)) {
	sw <- 1/sqrt(w)
	y <- sw * y
	X <- sw * X
	Z <- sw * Z
}

#  transform to independent observations
X <- as.matrix(X)
Z <- as.matrix(Z)
mx <- nrow(X)
nx <- ncol(X)
s <- La.svd(X,nu=mx,nv=0)
zeroeig <- abs(s$d/s$d[1]) < 1e-15
if(any(zeroeig))
	zero1 <- min((1:nx)[zeroeig])
else
	zero1 <- nx+1
Q <- s$u[,zero1:mx]
mq <- ncol(Q)
s <- La.svd(crossprod(Q,Z),nu=mq,nv=0)
uqy <- crossprod(s$u,(crossprod(Q,y)))
d <- rep(0,mq)
d[1:length(s$d)] <- s$d^2
dx <- cbind(Residual=1,Block=d)
dy <- uqy^2

#  low dimension cases
if(nrow(dx)==1) {
	sigma <- rep(NA,2)
	sigma[1] <- mean(dy)
	return(list(sigmasquared=sigma))
}
if(nrow(dx)==2) {
	sigma <- solve(dx,dy)
	v <- dx %*% sigma
	
}

#  fit gamma glm identity link to dy with dx as covariates
dfit <- glm.fit(dx,dy,family=Gamma(link="identity"),start=c(mean(dy),0))
dse <- sqrt(2*diag(chol2inv(dfit$qr$qr)))
out <- list(sigmasquared=dfit$coef,se.sigmasquared=dse)

#  fixed effect estimates
if(fixed.estimates) {
	v <- dfit$fitted.values
	s <- La.svd(Z,nu=mx,nv=0)
	d <- rep(0,mx)
	d[1:length(s$d)] <- s$d^2
	v <- drop( cbind(Residual=1,Block=d) %*% dfit$coefficients )
	mfit <- lm.wfit(crossprod(s$u,X),crossprod(s$u,y),1/v)
	out <- c(out,mfit)
	out$se.coefficients <- sqrt(diag(chol2inv(mfit$qr$qr)))
}

out
}
