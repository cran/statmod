\name{invgauss}
\alias{InverseGaussian}
\alias{dinvgauss}
\alias{pinvgauss}
\alias{qinvgauss}
\alias{rinvgauss}
\title{Inverse Gaussian Distribution}

\description{
Density, cumulative probability, quantiles and random generation for the inverse Gaussian distribution.
}

\usage{
dinvgauss(x, mean=1, shape=NULL, dispersion=1, log=FALSE)
pinvgauss(q, mean=1, shape=NULL, dispersion=1, lower.tail=TRUE, log.p=FALSE)
qinvgauss(p, mean=1, shape=NULL, dispersion=1, lower.tail=TRUE, log.p=FALSE,
          maxit=200L, tol=1e-14, trace=FALSE)
rinvgauss(n, mean=1, shape=NULL, dispersion=1)
}

\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{sample size. If \code{length(n)} is larger than 1, then \code{length(n)} random values are returned.}
  \item{mean}{vector of (positive) means.}
  \item{shape}{vector of (positive) shape parameters.}
  \item{dispersion}{vector of (positive) dispersion parameters. Ignored if \code{shape} is not \code{NULL}, in which case \code{dispersion=1/shape}.}
  \item{lower.tail}{logical; if \code{TRUE}, probabilities are P(X<q) otherwise P(X>q).}
  \item{log}{logical; if \code{TRUE}, the log-density is returned.}
  \item{log.p}{logical; if \code{TRUE}, probabilities are on the log-scale.}
  \item{maxit}{maximum number of Newton iterations used to find \code{q}.}
  \item{tol}{small positive numeric value giving the convergence tolerance for the quantile.}
  \item{trace}{logical, if \code{TRUE} then the working estimate for \code{q} from each iteration will be output.}
}

\value{
Output values give density (\code{dinvgauss}), probability (\code{pinvgauss}), quantile (\code{qinvgauss}) or random sample (\code{rinvgauss}) for the inverse Gaussian distribution with mean \code{mean} and dispersion \code{dispersion}.
Output is a vector of length equal to the maximum length of any of the arguments \code{x}, \code{q}, \code{mean}, \code{shape} or \code{dispersion}.
If the first argument is the longest, then all the attributes of the input argument are preserved on output, for example, a matrix \code{x} will give a matrix on output.
Elements of input vectors that are missing will cause the corresponding elements of the result to be missing, as will non-positive values for \code{mean} or \code{dispersion}.
}

\details{
The inverse Gaussian distribution takes values on the positive real line.
It is somewhat more right skew than the gamma distribution, with variance given by \code{dispersion*mean^3}.
The distribution has applications in reliability and survival analysis and is one of the response distributions used in generalized linear models.

Giner and Smyth (2016) show that the inverse Gaussian distribution converges to an inverse chi-squared distribution as the mean becomes large.

The functions provided here implement numeric algorithms developed by Giner and Smyth (2016) that achieve close to full machine accuracy for all possible parameter values.
Giner and Smyth (2016) show that the probability calculations provided by these functions are considerably more accurate, and in most cases faster, than previous implementations of inverse Gaussian functions.
The improvement in accuracy is most noticeable for extreme probability values and for large parameter values.

The shape and dispersion parameters are alternative parametrizations for the variability, with \code{dispersion=1/shape}.
Only one of these two arguments needs to be specified.
If both are set, then \code{shape} takes precedence.
}
 
\references{
Giner, G., and Smyth, G. K. (2016).
statmod: Probability calculations for the inverse Gaussian distribution.
\emph{R Journal} 8(1), 339-351.
\url{https://journal.r-project.org/archive/2016-1/giner-smyth.pdf}
}

\author{Gordon Smyth.

Very early S-Plus versions of these functions, using simpler algorithms, were published 1998 at \url{http://www.statsci.org/s/invgauss.html}.
Paul Bagshaw (Centre National d'Etudes des Telecommunications, France) contributed the original version of \code{qinvgauss} in December 1998.
Trevor Park (Department of Statistics, University of Florida) suggested improvements to a version of \code{rinvguass} in 2005.}

\examples{
q <- rinvgauss(10, mean=1, disp=0.5) # generate vector of 10 random numbers
p <- pinvgauss(q, mean=1, disp=0.5) # p should be uniformly distributed

# Quantile for small right tail probability:
qinvgauss(1e-20, mean=1.5, disp=0.7, lower.tail=FALSE)

# Same quantile, but represented in terms of left tail probability on log-scale
qinvgauss(-1e-20, mean=1.5, disp=0.7, lower.tail=TRUE, log.p=TRUE)
}

\keyword{distribution}
