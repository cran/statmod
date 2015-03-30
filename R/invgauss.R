dinvgauss <- function(x, mean=1, shape=NULL, dispersion=1, log=FALSE)
#	Probability density function of inverse Gaussian distribution
#	Gordon Smyth
#	Created 15 Jan 1998.  Last revised 30 March 2015.
{
#	Dispersion is reciprocal of shape
	if(!is.null(shape)) dispersion <- 1/shape

#	Make arguments same length
	r <- range(length(x),length(mean),length(dispersion))
	if(r[1L]==0L) return(numeric(0L))
	n <- r[2L]
	x <- rep_len(x,n)
	mu <- rep_len(mean,n)
	phi <- rep_len(dispersion,n)

#	Special cases
	NA.cases <- (is.na(x) | is.na(mu) | is.na(phi) | mu<=0 | phi<=0)
	left.limit <- x<=0
	right.limit <- x==Inf

#	Check for attributes
	has.attr <- !is.null(attributes(x))

	any.special <- has.attr | any(NA.cases) | any(left.limit) | any(right.limit)
	if(any.special) {
		logd <- x
		logd[left.limit] <- -Inf
		logd[right.limit] <- -Inf
		logd[NA.cases] <- NA
		ok <- !(NA.cases | left.limit | right.limit)
		logd[ok] <- .dinvgauss(x[ok],mean=mu[ok],dispersion=phi[ok],log=TRUE)
	} else {
		logd <- .dinvgauss(x,mean=mu,dispersion=phi,log=TRUE)
	}

	if(log) logd else exp(logd)
}

.dinvgauss <- function(x, mean=NULL, dispersion=1, log=FALSE)
#	Probability density function of inverse Gaussian distribution
#	with no argument checking and assuming mean=1
{
	notnullmean <- !is.null(mean)
	if(notnullmean) {
		x <- x/mean
		dispersion <- dispersion*mean
	}
	d <- (-log(dispersion)-log(2*pi)-3*log(x) - (x-1)^2/dispersion/x)/2
	if(notnullmean) d <- d-log(mean)
	if(log) d else exp(d)
}

pinvgauss <- function(q, mean=1, shape=NULL, dispersion=1, lower.tail=TRUE, log.p=FALSE)
#	Cumulative distribution function of inverse Gaussian distribution
#	Gordon Smyth
#	Created 15 Jan 1998.  Last revised 30 March 2015.
{
#	Dispersion is reciprocal of shape
	if(!is.null(shape)) dispersion <- 1/shape

#	Make arguments same length
	r <- range(length(q),length(mean),length(dispersion))
	if(r[1L]==0L) return(numeric(0L))
	n <- r[2L]
	q <- rep_len(q,n)
	mu <- rep_len(mean,n)
	phi <- rep_len(dispersion,n)

	nq <- length(q)
	if(nq==0) return(numeric(0))
	n <- max(nq,length(mean),length(dispersion))
	if(n>nq) q <- rep_len(q,n)
	mu <- rep_len(mean,n)
	phi <- rep_len(dispersion,n)

#	Special cases
	NA.cases <- (is.na(q) | is.na(mu) | is.na(phi) | mu<=0 | phi<=0)
	left.limit <- q<=0
	right.limit <- q==Inf

#	Check for attributes
	has.attr <- !is.null(attributes(q))

	any.special <- has.attr | any(NA.cases) | any(left.limit) | any(right.limit)
	if(any.special) {
		logp <- q
		if(lower.tail) {
			logp[left.limit] <- -Inf
			logp[right.limit] <- 0
		} else {
			logp[left.limit] <- 0
			logp[right.limit] <- -Inf
		}
		logp[NA.cases] <- NA
		ok <- !(NA.cases | left.limit | right.limit)
		logp[ok] <- .pinvgauss(q[ok],mean=mu[ok],dispersion=phi[ok],lower.tail=lower.tail,log.p=TRUE)
	} else {
		logp <- .pinvgauss(q,mean=mu,dispersion=phi,lower.tail=lower.tail,log.p=TRUE)
	}

	if(log.p) logp else(exp(logp))
}

.pinvgauss <- function(q, mean=NULL, dispersion=1, lower.tail=TRUE, log.p=FALSE)
#	Cumulative distribution function of inverse Gaussian distribution
#	without argument checking
{
	if(!is.null(mean)) {
		q <- q/mean
		dispersion <- dispersion*mean
	}
	pq <- sqrt(dispersion*q)
	a <- pnorm((q-1)/pq,lower.tail=lower.tail,log.p=TRUE)
	b <- 2/dispersion + pnorm(-(q+1)/pq,log.p=TRUE)
	if(lower.tail) b <- exp(b-a) else b <- -exp(b-a)
	logp <- a+log1p(b)
	if(log.p) logp else exp(logp)
}

rinvgauss <- function(n, mean=1, shape=NULL, dispersion=1)
#	Random variates from inverse Gaussian distribution
#	Gordon Smyth (with a correction by Trevor Park 14 June 2005)
#	Created 15 Jan 1998.  Last revised 30 March 2015.
{
#	Dispersion is reciprocal of shape
	if(!is.null(shape)) dispersion <- 1/shape

#	Check n
	if(length(n)>1L) n <- length(n) else n <- as.integer(n)
	if(n<0L) stop("n can't be negative")
	if(n==0L || length(mean)==0L || length(dispersion)==0L) return(numeric(0L))

#	Make arguments same length
	mu <- rep_len(mean,n)
	phi <- rep_len(dispersion,n)

#	Setup output vector
	r <- rep_len(0,n)

#	Non-positive parameters give NA
	i <- (mu > 0 & phi > 0)
	i[is.na(i)] <- FALSE
	if(!all(i)) {
		r[!i] <- NA
		n <- sum(i)
	}

#	Divide out mu	
	phi[i] <- phi[i]*mu[i]

	Y <- rchisq(n,df=1)
	X1 <- 1 + phi[i]/2 * (Y - sqrt(4*Y/phi[i]+Y^2))
	firstroot <- as.logical(rbinom(n,size=1L,prob=1/(1+X1)))
	r[i][firstroot] <- X1[firstroot]
	r[i][!firstroot] <- 1/X1[!firstroot]

	mu*r
}

qinvgauss  <- function(p, mean=1, shape=NULL, dispersion=1, lower.tail=TRUE, log.p=FALSE, maxit=40L, tol=1e-7, trace=FALSE)
#	Quantiles of the inverse Gaussian distribution
#
#	Current version using globally convergent Newton iteration
#	by created by Gordon Smyth 12 May 2014, last revised 30 March 2015.
#
#	Replaces an earlier function by Paul Bagshaw of 23 Dec 1998
{
#	Dispersion is reciprocal of shape
	if(!is.null(shape)) dispersion <- 1/shape

#	Make arguments same length
	r <- range(length(p),length(mean),length(dispersion))
	if(r[1L]==0L) return(numeric(0L))
	n <- r[2L]
	p <- rep_len(p,n)
	mu <- rep_len(mean,n)
	phi <- rep_len(dispersion,n)

#	Setup output
	q <- p

#	Special cases
	if(log.p) {
		NA.cases <- (is.na(p) | is.na(mu) | is.na(phi) | p>0 | mu<=0 | phi<=0)
		if(lower.tail) {
			left.limit <- p == -Inf
			right.limit <- p == 0
		} else {
			left.limit <- p == 0
			right.limit <- p == -Inf
		}
	} else {
		NA.cases <- (is.na(p) | is.na(mu) | is.na(phi) | p<0 | p>1 | mu<=0 | phi<=0)
		if(lower.tail) {
			left.limit <- p == 0
			right.limit <- p == 1
		} else {
			left.limit <- p == 1
			right.limit <- p == 0
		}
	}
	q[left.limit] <- 0
	q[right.limit] <- Inf
	q[NA.cases] <- NA
	ok <- !(NA.cases | left.limit | right.limit)

#	Convert to mean=1
	phi <- phi[ok]*mu[ok]
	p <- p[ok]

#	Mode of density and point of inflexion of cdf
	phi2 <- 1.5*phi
	x <- sqrt(1+phi2^2)-phi2
	if(trace) cat("mode",x,"\n")

	if(log.p) {
		step <- function(p,x,phi) {
			logF <- .pinvgauss(x,dispersion=phi,lower.tail=lower.tail,log.p=TRUE)
			logf <- .dinvgauss(x,dispersion=phi,log=TRUE)
			dlogp <- p-logF
			pos <- p>logF
			maxlogp <- logF
			maxlogp[pos] <- p[pos]
			d <- exp(maxlogp+log1p(-exp(-abs(dlogp)))-logf)
			d[!pos] <- -d[!pos]
			d
		}
	} else {
		step <- function(p,x,phi) (p - .pinvgauss(x, dispersion=phi, lower.tail=lower.tail)) / .dinvgauss(x, dispersion=phi)
	}

#	Newton iteration is monotonically convergent from point of inflexion
	iter <- 0
	i <- rep_len(TRUE,length(phi))
	while(any(i)) {
		iter <- iter+1
		if(iter > maxit) {
			warning("max iterations exceeded")
			break
		}
		dx <- step(p[i],x[i],phi[i])
		if(lower.tail)
			x[i] <- x[i] + dx
		else
			x[i] <- x[i] - dx
		i[i] <- (abs(dx) > tol)
		if(trace) {
			cat("Iter=",iter,"Still converging=",sum(i),"\n")
			if(n < 6L)
				cat("x ",x,"\ndx ",dx,"\n")
			else
				cat("Quantiles x ",quantile(x),"\nMax dx ",max(dx),"\n")
		}
	}

#	Mu scales the distribution
	q[ok] <- x*mu[ok]
	q
}
