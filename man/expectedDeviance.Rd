\name{expectedDeviance}
\alias{expectedDeviance}
\title{Expected Value of Scaled Unit Deviance for Linear Exponential Families}

\description{
Expected value and variance of the scaled unit deviance for common generalized linear model families.
}

\usage{
expectedDeviance(mu, family="binomial", binom.size, nbinom.size, gamma.shape)
}

\arguments{
  \item{mu}{numeric vector or matrix giving mean of response variable.}
  \item{family}{character string indicating the linear exponential family. Possible values are \code{"binomial"},\code{"gaussian"}, \code{"Gamma"}, \code{"inverse.gaussian"}, \code{"poisson"} or \code{"negative.binomial"}.}
  \item{binom.size}{integer vector giving the number of binomial trials when \code{family = "binomial"}. Equivalent to the \code{"size"} argument of \code{pbinom}.}
  \item{nbinom.size}{numeric vector giving the negative binomial size parameter when \code{family = "negative.binomial"}, such that the variance of the response variable is \code{mu + mu^2 / nbinom.size}. Equivalent to the \code{"size"} parameter of \code{pnbinom}.}
  \item{gamma.shape}{numeric vector giving the gamma shape parameter when \code{family = "Gamma"}, such that the variance of the response variable is \code{mu^2 / gamma.shape}. Equivalent to the \code{"shape"} parameter of \code{pgamma}.}
}

\details{
For a generalized linear model (GLM), the scaled unit deviances can be computed using \code{d <- f$dev.resids(y, mu, wt=1/phi)} where \code{f} is the GLM family object, \code{y} is the response variable, \code{mu} is the vector of means and \code{phi} is the vector of GLM dispersions (incorporating any prior weights).

The scaled unit deviances are often treated as being chiquare distributed on 1 df, so the mean should be 1 and the variance should be 2.
This distribution result only holds however when the saddlepoint approximation is accurate for the response variable distribution (Dunn and Smyth, 2018).
In other cases, the expected value and variance of the unit deviances can be far from the nominal values.
The \code{expectedDeviance} function returns the exact mean and variance of the unit deviance for the usual GLM familes assuming that \code{mu} is the true mean and \code{phi} is the true dispersion.

When \code{family} is \code{"poisson"}, \code{"binomial"} or \code{"negative.binomial"}, the expected values and variances are computed using Chebyshev polynomial approximations.
When \code{family = "Gamma"}, the function uses exact formulas derived by Smyth (1989).
}

\value{
A list with the components
  \item{mean}{expected values}
  \item{variance}{variances}
both of which have the same length and dimensions as the input \code{mu}.
}

\author{Lizong Chen and Gordon Smyth}

\references{
Dunn PK, Smyth GK (2018).
\emph{Generalized linear models with examples in R}.
Springer, New York, NY.
\doi{10.1007/978-1-4419-0118-7}

Smyth, G. K. (1989).
Generalized linear models with varying dispersion.
\emph{J. R. Statist. Soc. B}, \bold{51}, 47-61.
\doi{10.1111/j.2517-6161.1989.tb01747.x}
}

\examples{
# Poisson example
lambda <- 3
nsim <- 1e4
y <- rpois(nsim, lambda=lambda)
d <- poisson()$dev.resids(y=y, mu=rep(lambda,nsim), wt=1)
c(mean=mean(d), variance=var(d))
unlist(expectedDeviance(mu=lambda, family="poisson"))

# binomial example
n <- 10
p <- 0.01
y <- rbinom(nsim, prob=p, size=n)
d <- binomial()$dev.resids(y=y/n, mu=rep(p,nsim), wt=n)
c(mean=mean(d), variance=var(d))
unlist(expectedDeviance(mu=p, family="binomial", binom.size=n))

# gamma example
alpha <- 5
beta <- 2
y <- beta * rgamma(1e4, shape=alpha)
d <- Gamma()$dev.resids(y=y, mu=rep(alpha*beta,n), wt=alpha)
c(mean=mean(d), variance=var(d))
unlist(expectedDeviance(mu=alpha*beta, family="Gamma", gamma.shape=alpha))

# negative binomial example
library(MASS)
mu <- 10
phi <- 0.2
y <- rnbinom(nsim, mu=mu, size=1/phi)
f <- MASS::negative.binomial(theta=1/phi)
d <- f$dev.resids(y=y, mu=rep(mu,nsim), wt=1)
c(mean=mean(d), variance=var(d))
unlist(expectedDeviance(mu=mu, family="negative.binomial", nbinom.size=1/phi))

# binomial expected deviance tends to zero for p small:
p <- seq(from=0.001,to=0.11,len=200)
ed <- expectedDeviance(mu=p,family="binomial",binom.size=10)
plot(p,ed$mean,type="l")
}

\seealso{
\code{\link{family}}, \code{\link{meanval.digamma}}, \code{\link{d2cumulant.digamma}}.
}

\keyword{distributions}
