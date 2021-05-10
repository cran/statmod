\name{tweedie}
\alias{tweedie}
\title{Tweedie Generalized Linear Models}

\description{
Produces a generalized linear model family object with any power variance function and any power link. Includes the Gaussian, Poisson, gamma and inverse-Gaussian families as special cases.
}
\usage{
tweedie(var.power = 0, link.power = 1 - var.power)
}

\arguments{
\item{var.power}{index of power variance function}
\item{link.power}{index of power link function. \code{link.power=0} produces a log-link. Defaults to the canonical link, which is \code{1-var.power}.} 
}

\value{
A family object, which is a list of functions and expressions used by \code{glm} and \code{gam} in their iteratively reweighted least-squares algorithms.
See \code{\link{family}} and \code{\link{glm}} in the R base help for details. 
}

\details{
This function provides access to a range of generalized linear model (GLM) response distributions that are not otherwise provided by R.
It is also useful for accessing distribution/link combinations that are disallowed by the R \code{glm} function.
The variance function for the GLM is assumed to be V(mu) = mu^var.power, where mu is the expected value of the distribution.
The link function of the GLM is assumed to be mu^link.power for non-zero values of link.power or log(mu) for var.power=0.
For example, \code{var.power=1} produces the identity link.
The canonical link for each Tweedie family is \code{link.power = 1 - var.power}.

The Tweedie family of GLMs is discussed in detail by Dunn and Smyth (2018).
Each value of \code{var.power} corresponds to a particular type of response distribution.
The values 0, 1, 2 and 3 correspond to the normal distribution, the Poisson distribution, the gamma distribution and the inverse-Gaussian distribution respectively.
For these choices of \code{var.power}, the Tweedie family is exactly equivalent to the usual GLM famly except with a greater choice of link powers.
For example, \code{tweedie(var.power = 1, link.power = 0)} is exactly equivalent to \code{poisson(link = "log")}.

The most interesting Tweedie families occur for \code{var.power} between 1 and 2.
For these GLMs, the response distribution has mass at zero (i.e., it has exact zeros) but is otherwise continuous on the positive real numbers (Smyth, 1996; Hasan et al, 2012).
These GLMs have been used to model rainfall for example.
Many days there is no rain at all (exact zero) but, if there is any rain, then the actual amount of rain is continuous and positive.

Generally speaking, \code{var.power} should be chosen so that the theoretical response distribution matches the type of response data being modeled.
Hence \code{var.power} should be chosen between 1 and 2 only if the response observations are continuous and positive except for exact zeros and \code{var.power} should be chosen greater than or equal to 2 only if the response observations are continuous and strictly positive.

There are no theoretical Tweedie GLMs with var.power between 0 and 1 (Jorgensen 1987).
The \code{tweedie} function will work for those values but the family should be interpreted in a quasi-likelihood sense.

Theoretical Tweedie GLMs do exist for negative values of var.power, but they are of little practical application.
These distributions assume 
The \code{tweedie} function will work for those values but the family should be interpreted in a quasi-likelihood sense.

The name Tweedie has been associated with this family by Joergensen (1987) in honour of M. C. K. Tweedie.
Joergensen (1987) gives a mathematical derivation of the Tweedie distributions proving that no distributions exist for var.power between 0 and 1.

Mathematically, a Tweedie GLM assumes the following.
Let \eqn{\mu_i = E(y_i)} be the expectation of the \eqn{i}th response. We assume that
\deqn{\mu_i^q = x_i^Tb, var(y_i) = \phi \mu_i^p} 

where \eqn{x_i} is a vector of covariates and b is a vector of regression cofficients, for some \eqn{\phi}, \eqn{p} and \eqn{q}.
This family is specified by \code{var.power = p} and \code{link.power = q}.
A value of zero for \eqn{q} is interpreted as \eqn{\log(\mu_i) = x_i^Tb}. 

The following table summarizes the possible Tweedie response distributions:

\tabular{cl}{
\bold{var.power} \tab \bold{Response distribution}\cr
0 \tab Normal\cr
1 \tab Poisson\cr
(1, 2) \tab Compound Poisson, non-negative with mass at zero\cr
2 \tab Gamma\cr
3 \tab Inverse-Gaussian\cr
> 2 \tab Stable, with support on the positive reals
}

}
 
\references{
Dunn, P. K., and Smyth, G. K, (2018).
\emph{Generalized linear models with examples in R}. Springer, New York, NY.
\doi{10.1007/978-1-4419-0118-7}
(Chapter 12 gives an overall discussion of Tweedie GLMs with R code and case studies.)

Hasan, M.M. and Dunn, P.K. (2012).
Understanding the effect of climatology on monthly rainfall amounts in Australia using Tweedie GLMs.
\emph{International Journal of Climatology}, 32(7) 1006-1017.
(An example with var.power between 1 and 2)

Joergensen, B. (1987). Exponential dispersion models. \emph{J. R. Statist. Soc.} B \bold{49}, 127-162. 
(Mathematical derivation of Tweedie response distributions)

Tweedie, M. C. K. (1984). An index which distinguishes between some important exponential families.
In \emph{Statistics: Applications and New Directions}. Proceedings of the Indian Statistical Institute Golden Jubilee International Conference.
(Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
(The original mathematical paper from which the family is named)

Smyth, G. K. (1996). Regression modelling of quantity data with exact zeroes.
\emph{Proceedings of the Second Australia-Japan Workshop on Stochastic Models in Engineering, Technology and Management}.
Technology Management Centre, University of Queensland, pp. 572-580.
\url{http://www.statsci.org/smyth/pubs/RegressionWithExactZerosPreprint.pdf}
(Derivation and examples of Tweedie GLMS with var.power between 0 and 1)

Smyth, G. K., and Verbyla, A. P., (1999). Adjusted likelihood methods for modelling dispersion in generalized linear models. \emph{Environmetrics} \bold{10}, 695-709.
\url{http://www.statsci.org/smyth/pubs/Ties98-Preprint.pdf}
(Includes examples of Tweedie GLMs with \code{var.power=2} and \code{var.power=4})
}

\author{Gordon Smyth}

\seealso{\code{\link{glm}}, \code{\link{family}}, \code{\link[tweedie]{dtweedie}}}

\examples{
y <- rgamma(20,shape=5)
x <- 1:20
# Fit a poisson generalized linear model with identity link
glm(y~x,family=tweedie(var.power=1,link.power=1))

# Fit an inverse-Gaussion glm with log-link
glm(y~x,family=tweedie(var.power=3,link.power=0)) 
}

\keyword{regression}
