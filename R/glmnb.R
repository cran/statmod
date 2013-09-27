glmnb.fit <- function(X,y,dispersion,offset=0,coef.start=NULL,start.method="mean",tol=1e-6,maxit=50,trace=FALSE)
#  Fit negative binomial generalized linear model with log link
#  by Levenberg damped Fisher scoring
#  Yunshun Chen and Gordon Smyth
#  2 November 2010.  Last modified 20 November 2012.
{
#  Check input values for y
	y <- as.vector(y)
	if(any(y < 0)) stop("y must be non-negative")
	if(!all(is.finite(y))) stop("All y values must be finite and non-missing")
	ymax <- max(y)
	n <- length(y)

#	Handle zero length y as special case
	if(n == 0) stop("y has length zero")

#	Check input values for X
	X <- as.matrix(X)
	if(n != nrow(X)) stop("length(y) not equal to nrow(X)")
	if(!all(is.finite(X))) stop("All X values must be finite and non-missing")
	p <- ncol(X)
	if(p > n) stop("More columns than rows in X")

#	Handle y all zero as special case
	if(ymax==0) return(list(coefficients=rep(0,p),fitted.values=rep(0,n),deviance=0,iter=0,convergence="converged"))

#	Check input values for dispersion
	if(any(dispersion<0)) stop("dispersion values must be non-zero")
	phi <- rep(dispersion,length.out=n)

#	Check input values for offset
	offset <- rep(offset,length=n)

#  Starting values
	delta <- min(ymax,1/6)
	y1 <- pmax(y,delta)
	if(is.null(coef.start)) {
		start.method <- match.arg(start.method,c("log(y)","mean"))
		if(start.method=="log(y)") {
			fit <- lm.fit(X,log(y1)-offset)
			beta <- fit$coefficients
			mu <- exp(fit$fitted.values+offset)
		} else {
			N <- exp(offset)
			rate <- y/N
			w <- N/(1+phi*N)
			beta.mean <- log(sum(w*rate)/sum(w))
			beta <- qr.coef(qr(X),rep(beta.mean,length=n))
			mu <- drop(exp(X %*% beta + offset))
		}
	} else {
		beta <- coef.start
		mu <- drop(exp(X %*% beta + offset))
	}

	unit.dev.poisson <- function(y,mu) {
		2 * ( y*log(y/mu)- (y-mu) )
	}
	unit.dev.gamma <- function(y,mu) {
		2 * ( (y-mu)/mu - log(y/mu))
	}
	unit.dev.negbin <- function(y,mu,size) {
		2 * ( y*log(y/mu) - (y+size)*log((y+size)/(mu+size)) )
	}
	total.deviance <- function(y,mu,phi) {
		if(any(is.infinite(mu))) return(Inf)
		poisson.like <- (phi*mu < 1e-6)
		gamma.like <- (phi*mu > 1e6)
		negbin <- !(poisson.like | gamma.like)
		y <- y+1e-8
		mu <- mu+1e-8
		unit.dev <- y
		if(any(poisson.like)) unit.dev[poisson.like] <- unit.dev.poisson(y[poisson.like],mu[poisson.like])
		if(any(gamma.like)) {
			m <- mu[gamma.like]
			alpha <- m/(1+phi[gamma.like]*m)
			unit.dev[gamma.like] <- unit.dev.gamma(y[gamma.like],m)*alpha
		}
		if(any(negbin)) unit.dev[negbin] <- unit.dev.negbin(y[negbin],mu[negbin],1/phi[negbin])
		sum(unit.dev)
	}

	dev <- total.deviance(y,mu,phi)

#	Scoring iteration with Levenberg damping
	iter <- 0
	if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")
	repeat {
		iter <- iter+1

#		test for iteration limit
		if(iter > maxit) break

#		information matrix
		v.div.mu <- 1+phi*mu
		XVX <- crossprod(X,(mu/v.div.mu)*X)
		maxinfo <- max(diag(XVX))
		if(iter==1) {
			lambda <- maxinfo * 1e-6
			lambda <- max(lambda,1e-13)
			lambdaceiling <- maxinfo * 1e13
			lambdabig <- FALSE
			I <- diag(p)
		}

#		score vector
		dl <- crossprod(X,(y-mu)/v.div.mu)

#		Levenberg damping
		dbeta <- beta
		lev <- 0
		repeat {
			lev <- lev+1

#			trial step
			R <- chol(XVX + lambda*I, pivot=TRUE)
			while(attr(R,"rank")<p) {
				lambda <- lambda*10
				R <- chol(XVX + lambda*I, pivot=TRUE)
			}
			j <- attr(R,"pivot")
			dbeta[j] <- backsolve(R,backsolve(R,dl[j],transpose=TRUE))
			betanew <- beta + dbeta
			munew <- drop(exp(X %*% betanew + offset))
			devnew <- total.deviance(y,munew,phi)
			if(devnew <= dev || devnew/ymax < 1e-13) {
				beta <- betanew
				mu <- munew
				dev <- devnew
				break
			}

#			step not successful so increase damping
			lambda <- 2*lambda
			if(trace) cat("Damping increased to",lambda,"\n")

#			exit if too much damping
			lambdabig <- (lambda > lambdaceiling)
			if(lambdabig) {
				warning("Too much damping - convergence tolerance not achievable")
				break
			}
		}

#		iteration output
		if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")

#		keep exiting if too much damping
		if(lambdabig) break

#		test for convergence
		scoresquare <- crossprod(dl,dbeta)
		if(trace) cat("Convergence criterion",scoresquare,dl,dbeta,"\n")
		if( scoresquare < tol || dev/ymax < 1e-12) break

#		decrease damping if successful at first try
		if(lev==1) lambda <- lambda/10
	}

	beta <- drop(beta)
	names(beta) <- colnames(X)
	convergence <- "converged"
	if(lambdabig) convergence <- "lambdabig"
	if(iter>maxit) convergence <- "maxit"
	list(coefficients=beta,fitted.values=as.vector(mu),deviance=dev,iter=iter,convergence=convergence)
}

