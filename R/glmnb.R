####################################
########### glmnb.fit ##############
####################################

glmnb.fit <- function(X,y,dispersion,offset=0,start=NULL,tol=1e-6,maxit=50,trace=FALSE)
#  Fit negative binomial generalized linear model with log link
#  by Levenberg damped Fisher scoring
#  Yunshun Chen and Gordon Smyth
#  2 November 2010.  Last modified 3 November 2010.
{
#  check input
	X <- as.matrix(X)
	n <- nrow(X)
	p <- ncol(X)
	if(p > n) stop("More columns than rows in X")
	y <- as.vector(y)
	if(n != length(y)) stop("length(y) not equal to nrow(X)")
	if(n == 0) return(list(coefficients=numeric(0),fitted.values=numeric(0),deviance=numeric(0)))
	if(!(all(is.finite(y)) || all(is.finite(X)))) stop("All values must be finite and non-missing")
	if(any(y < 0)) stop("y must be non-negative")
	maxy <- max(y)
	if(maxy==0) return(list(coefficients=rep(0,p),fitted.values=rep(0,n),deviance=NA))
	y1 <- pmax(y,1/6)
	phi <- dispersion

#  starting values
	if(is.null(start)) {
		fit <- lm.fit(X,log(y1)-offset)
		beta <- fit$coefficients
		mu <- exp(fit$fitted.values+offset)
	} else {
		beta <- start
		mu <- exp(X %*% beta + offset)
	}

	deviance.nb <- function(y,mu,phi) {
		if(any(mu<0)) return(Inf)
		o <- (y < 1e-14) & (mu < 1e-14)
		if(any(o)) {
			if(all(o)) {
				dev <- 0
			} else {
				y1 <- y[!o]
				mu1 <- mu[!o]
				dev <- 2*sum(y1*log(y1/mu1) + (y1+1/phi)*log((mu1+1/phi)/(y1+1/phi)) )
			}
		} else {
			dev <- 2*sum(y*log(y1/mu) + (y+1/phi)*log((mu+1/phi)/(y+1/phi)) )
		}
	}

	dev <- deviance.nb(y,mu,phi)

#	Scoring iteration with Levenberg damping
	iter <- 0
	if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")
	repeat {
		iter <- iter+1

#		information matrix
		v.div.mu <- 1+phi*mu
		XVX <- crossprod(X,(mu/v.div.mu)*X)
		maxinfo <- max(diag(XVX))
		if(iter==1) {
			lambda <- abs(mean(diag(XVX)))/p
			I <- diag(p)
		}

#		score vector
		dl <- crossprod(X,(y-mu)/v.div.mu)

#		Levenberg damping
		betaold <- beta
		devold <- dev
		lev <- 0
		repeat {
			lev <- lev+1

#			trial step
			R <- chol(XVX + lambda*I)
			dbeta <- backsolve(R,backsolve(R,dl,transpose=TRUE))
			beta <- betaold + dbeta
			mu <- exp(X %*% beta + offset)
			dev <- deviance.nb(y,mu,phi)
			if(dev <= devold || dev/max(mu) < 1e-13) break

#			exit if too much damping
			if(lambda/maxinfo > 1e13) {
				beta <- betaold
				warning("Too much damping - convergence tolerance not achievable")
				break
			}

#			step not successful so increase damping
			lambda <- 2*lambda
			if(trace) cat("Damping increased to",lambda,"\n")
		}

#		iteration output
		if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")

#		keep exiting if too much damping
		if(lambda/maxinfo > 1e13) break

#		decrease damping if successful at first try
		if(lev==1) lambda <- lambda/10

#		test for convergence
		if( crossprod(dl,dbeta) < tol || dev/max(mu) < 1e-13) break

#		test for iteration limit
		if(iter > maxit) break
	}

	beta <- drop(beta)
	names(beta) <- colnames(X)
	list(coefficients=beta,fitted.values=as.vector(mu),deviance=dev,iter=iter)
}

