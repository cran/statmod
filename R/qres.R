qres.binom <- function(glm.obj)
{
# Randomized quantile residuals for binomial glm
# GKS  20 Oct 96. Last modified 25 Jan 02.
#
	p <- fitted(glm.obj)
	y <- glm.obj$y
	if(!is.null(glm.obj$prior.weights))
		n <- glm.obj$prior.weights
	else
		n <- rep(1,length(y))
	y <- n * y
	a <- pbinom(y - 1, n, p)
	b <- pbinom(y, n, p)
	u <- runif(n = length(y), min = a, max = b)
	qnorm(u)
}

qres.pois <- function(glm.obj)
{
#	Quantile residuals for Poisson glm
#	GKS  28 Dec 96
#
	y <- glm.obj$y
	mu <- fitted(glm.obj)
	a <- ppois(y - 1, mu)
	b <- ppois(y, mu)
	u <- runif(n = length(y), min = a, max = b)
	qnorm(u)
}

qres.gamma <- function(glm.obj, dispersion = NULL)
{
# Quantile residuals for gamma glm
# GKS  28 Dec 96, 10 Jan 97
#
	mu <- fitted(glm.obj)
	y <- glm.obj$y
	df <- glm.obj$df.residual
	w <- glm.obj$prior.weights
	if(is.null(w))
		w <- 1
	if(is.null(dispersion))
		dispersion <- sum(w * ((y - mu)/mu)^2)/df
	u <- pgamma((w * y)/mu/dispersion, w/dispersion)
	qnorm(u)
}

qres.invgauss <- function(glm.obj, dispersion = NULL)
{
# Quantile residuals for inverse Gaussian glm
# GKS  15 Jan 98
#
	mu <- fitted(glm.obj)
	y <- glm.obj$y
	df <- glm.obj$df.residual
	w <- glm.obj$prior.weights
	if(is.null(w))
		w <- 1
	if(is.null(dispersion))
		dispersion <- sum(w * (y - mu)^2 / (mu^2*y)) / df
	u <- pinvgauss(y, mu, lambda=1/dispersion)
	qnorm(u)
}

qres.nbinom <- function(glm.obj)
{
#	Quantile residuals for Negative Binomial glm
#	GKS  22 Jun 97
#
	y <- glm.obj$y
	size <- glm.obj$call$family$a
	mu <- fitted(glm.obj)
	p <- size/(mu + size)
	a <- ifelse(y > 0, pbeta(p, size, pmax(y, 1)), 0)
	b <- pbeta(p, size, y + 1)
	u <- runif(n = length(y), min = a, max = b)
	qnorm(u)
}

qres.tweedie <- function(glm.obj, dispersion = NULL)
{
# Quantile residuals for Tweedie glms
# GKS  29 April 98.  Revised 20 April 99, 20 July 2001
#
	mu <- fitted(glm.obj)
	y <- glm.obj$y
	df <- glm.obj$df.residual
	w <- glm.obj$prior.weights
	if(is.null(w))
		w <- 1
	if(is.null(glm.obj$call$var.power))
		p <- 1.5
	else
		p <- eval(glm.obj$call$var.power, local = sys.parent())
	if(is.null(dispersion))
		dispersion <- sum((w * (y - mu)^2)/mu^p)/df
	u <- ptweedie(y, fitted(glm.obj), dispersion/w, p)
	if(p>1&&p<2)
		u[y == 0] <- runif(sum(y == 0), min = 0, max = u[y == 0])
	qnorm(u)
}
