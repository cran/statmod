\name{gauss.quad.prob}
\alias{gauss.quad.prob}
\title{Gaussian Quadrature with Probability Distributions}

\description{Calculate nodes and weights for Gaussian quadrature in terms of probability distributions.}

\usage{gauss.quad.prob(n,dist="uniform",l=0,u=1,mu=0,sigma=1,alpha=1,beta=1)}

\arguments{
  \item{n}{number of nodes and weights}
  \item{dist}{distribution that Gaussian quadrature is based on, one of \code{"uniform"}, \code{"normal"}, \code{"beta"} or \code{"gamma"}}
  \item{l}{lower limit of uniform distribution}
  \item{u}{upper limit of uniform distribution}
  \item{mu}{mean of normal distribution}
  \item{sigma}{standard deviation of normal distribution}
  \item{alpha}{positive shape parameter for gamma distribution or first shape parameter for beta distribution}
  \item{beta}{positive scale parameter for gamma distribution or second shape parameter for beta distribution}
}

\value{A list containing the components
  \item{nodes}{vector of values at which to evaluate the function}
  \item{weights}{vector of weights to give the function values}
}

\details{
This is a rewriting and simplification of \code{gauss.quad} in terms of probability distributions.

The expected value of \code{f(X)} is approximated by \code{sum(w*f(x))} where \code{x} is the vector of nodes and \code{w} is the vector of weights.  The approximation is exact if \code{f(x)} is a polynomial of order no more than \code{2n-1}.
The possible choices for the distribution of \code{X} are as follows:

Uniform on \code{(l,u)}.

Normal with mean \code{mu} and standard deviation \code{sigma}.

Beta with density \code{x^(alpha-1)*(1-x)^(beta-1)/B(alpha,beta)} on \code{(0,1)}.

Gamma with density \code{x^(alpha-1)*exp(-x/beta)/beta^alpha/gamma(alpha)}.
}

\references{
Golub, G. H., and Welsch, J. H. (1969). Calculation of Gaussian
quadrature rules. \emph{Mathematics of Computation} \bold{23}, 221-230.

Golub, G. H. (1973). Some modified matrix eigenvalue problems.
\emph{Siam Review} \bold{15}, 318-334.

Smyth, G. K. (1998). Polynomial approximation. In: \emph{Encyclopedia of Biostatistics}, P. Armitage and T. Colton (eds.), Wiley, London, pp. 3425-3429.
\url{http://www.statsci.org/smyth/pubs/poly.ps}

Stroud and Secrest (1966). \emph{Gaussian Quadrature Formulas}. Prentice-
Hall, Englewood Cliffs, N.J.
}

\author{Gordon Smyth}

\seealso{
\code{\link[statmod:gauss.quad]{gauss.quad}}, \code{\link{integrate}}
}

\examples{
out <- gauss.quad.prob(10,"normal")
sum(out$weights * out$nodes^4)
#  the 4th moment of the standard normal is 3

out <- gauss.quad.prob(32,"gamma",alpha=5)
sum(out$weights * log(out$nodes))
#  the expected value of log(X) where X is gamma is digamma(alpha)
}
\keyword{math}
