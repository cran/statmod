dinvgauss <- function(x, mu, lambda = 1, log=FALSE)
#  Density of inverse Gaussian distribution
#  Gordon Smyth
#	15 Jan 1998.  Last revised 19 June 2009.
{
	if(any(mu<=0)) stop("mu must be positive")
	if(any(lambda<=0)) stop("lambda must be positive")
	d <- rep(-Inf,length(x))
	i <- x>0
	d[i] <- (log(lambda)-log(2*pi)-3*log(x[i]))/2-lambda*(x[i]-mu)^2/(2*mu^2*x[i])
	if(!log) d <- exp(d)
	if(!is.null(Names <- names(x))) names(d) <- rep(Names, length = length(d))
	d
}

pinvgauss <- function(q, mu, lambda = 1)
{
#  Inverse Gaussian distribution function
#  GKS  15 Jan 98
#
	if(any(mu<=0)) stop("mu must be positive")
	if(any(lambda<=0)) stop("lambda must be positive")
	n <- length(q)
	if(length(mu)>1 && length(mu)!=n) mu <- rep(mu,length=n)
	if(length(lambda)>1 && length(lambda)!=n) lambda <- rep(lambda,length=n)
	lq <- sqrt(lambda/q)
	qm <- q/mu
	p <- ifelse(q>0,pnorm(lq*(qm-1))+exp(2*lambda/mu)*pnorm(-lq*(qm+1)),0)
	if(!is.null(Names <- names(q)))
		names(p) <- rep(Names, length = length(p))
	p
}

rinvgauss <- function(n, mu, lambda = 1)
#	Random variates from inverse Gaussian distribution
#	Reference:  Chhikara and Folks, The Inverse Gaussian
#	Distribution, Marcel Dekker, 1989, page 53.
#	Gordon Smyth 15 Jan 98.
#	Revised by Trevor Park 14 June 2005.
{
	if(any(mu<=0)) stop("mu must be positive")
	if(any(lambda<=0)) stop("lambda must be positive")
	if(length(n)>1) n <- length(n)
	if(length(mu)>1 && length(mu)!=n) mu <- rep(mu,length=n)
	if(length(lambda)>1 && length(lambda)!=n) lambda <- rep(lambda,length=n)
	y2 <- rchisq(n,1)
	u <- runif(n)
	r2 <- mu/(2*lambda)*(2*lambda+mu*y2+sqrt(4*lambda*mu*y2+mu^2*y2^2))
	r1 <- mu^2/r2
	ifelse(u < mu/(mu+r1), r1, r2)
}

qinvgauss  <- function(p, mu, lambda = 1, iter=20L, trace=FALSE)
#	Quantiles of the inverse Gaussian distribution

#	Original version by Dr Paul Bagshaw
#	Centre National d'Etudes des Telecommunications (DIH/DIPS)
#	Technopole Anticipa, France
#	paul.bagshaw@cnet.francetelecom.fr
#	23 Dec 1998

#	Current version by Gordon Smyth
#	13 April 2014
{
	if(any(mu <= 0)) stop("mu must be positive")
	if(any(lambda <= 0)) stop("lambda must be positive")
	n <- length(p)
	if(length(mu) > 1 && length(mu) != n) mu <- rep(mu, length = n)
	if(length(lambda) > 1 && length(lambda) != n) lambda <- rep(lambda, length = n)

#	Shape of distribution depends only on phi
	phi <- lambda / mu

#	Mode of density and point of inflexion of cdf (when mu=1)
	x <- sqrt(1+9/4/phi^2)-3/2/phi
	if(trace) cat("mode",x,"\n")

#	Newton iteration is monotonically convergent from point of inflexion
	for (i in 1L:iter) {
		cum <- pinvgauss (x, 1, phi)
		dx <- (cum - p) / dinvgauss (x, 1, phi)
		if (all(abs(dx) < 1e-5)) break
		x <- x - dx
		if(trace) cat(iter,x,"\n")
	}

#	Mu scales the distribution
	x * mu
}
