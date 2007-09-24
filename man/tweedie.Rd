\name{tweedie}
\alias{tweedie}
\title{Tweedie Generalized Linear Models}

\description{
Produces a generalized linear model family object with any power variance function and any power link. Includes the Gaussian, Poisson, gamma and inverse-Gaussian families as special cases.
}
\usage{
tweedie(var.power=0, link.power=1-var.power)
}

\arguments{
\item{var.power}{index of power variance function}
\item{link.power}{index of power link function. \code{link.power=0} produces a log-link. Defaults to the canonical link, which is \code{1-var.power}.} 
}

\value{
A family object, which is a list of functions and expressions used by glm and gam in their iteratively reweighted least-squares algorithms.
See \code{family} and \code{glm} in the R base help for details. 
}

\details{
This function provides access to a range of generalized linear model response distributions which are not otherwise provided by R, or any other package for that matter. It is also useful for accessing distribution/link combinations which are disallowed by the R \code{glm} function.

Let \eqn{\mu_i = E(y_i)} be the expectation of the \eqn{i}th response. We assume that
\deqn{\mu_i^q = x_i^Tb, var(y_i) = \phi \mu_i^p} 

where \eqn{x_i} is a vector of covariates and b is a vector of regression cofficients, for some \eqn{\phi}, \eqn{p} and \eqn{q}. This family is specified by \code{var.power = p} and \code{link.power = q}. A value of zero for \eqn{q} is interpreted as \eqn{\log(\mu_i) = x_i^Tb}. 

The variance power \eqn{p} characterizes the distribution of the responses \eqn{y}. The following are some special cases:

\tabular{cl}{
\bold{p} \tab \bold{Response distribution}\cr
0 \tab Normal\cr
1 \tab Poisson\cr
(1, 2) \tab Compound Poisson, non-negative with mass at zero\cr
2 \tab Gamma\cr
3 \tab Inverse-Gaussian\cr
> 2 \tab Stable, with support on the positive reals
}

The name Tweedie has been associated with this family by Joergensen (1987) in honour of M. C. K. Tweedie. 
}
 
\references{
Tweedie, M. C. K. (1984). An index which distinguishes between some important exponential families. In \emph{Statistics: Applications and New Directions}. Proceedings of the Indian Statistical Institute Golden Jubilee International Conference. (Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.

Joergensen, B. (1987). Exponential dispersion models. \emph{J. R. Statist. Soc.} B \bold{49}, 127-162. 

Smyth, G. K. (1996). Regression modelling of quantity data with exact zeroes. Proceedings of the Second Australia-Japan Workshop on Stochastic Models in Engineering, Technology and Management. Technology Management Centre, University of Queensland, pp. 572-580. 

Joergensen, B. (1997). \emph{Theory of Dispersion Models}, Chapman and Hall, London. 

Smyth, G. K., and Verbyla, A. P., (1999). Adjusted likelihood methods for modelling dispersion in generalized linear models. \emph{Environmetrics} \bold{10}, 695-709. 
}

\author{Gordon Smyth}

\seealso{\code{\link[base]{glm}}, \code{\link[base]{family}}, \code{\link[tweedie]{dtweedie}}}

\examples{
y <- rgamma(20,shape=5)
x <- 1:20
# Fit a poisson generalized linear model with identity link
glm(y~x,family=tweedie(var.power=1,link.power=1))

# Fit an inverse-Gaussion glm with log-link
glm(y~x,family=tweedie(var.power=3,link.power=0)) 
}

\keyword{regression}
