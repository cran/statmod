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
          maxit=40L, tol=1e-7, trace=FALSE)
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
The inverse Gaussian distribution takes values on the positive real line (Tweedie, 1957; Chhikara and Folks, 1989).
It is somewhat more right skew than the gamma distribution, with variance given by \code{dispersion*mean^3}.
The distribution has applications in reliability and survival analysis, and is one of the response distributions used in generalized linear models.

The shape and dispersion parameters are alternative parametrizations for the variability, with \code{dispersion=1/shape}.
Only one of these two arguments needs to be specified.
If both are set, then \code{shape} takes precedence.

\code{pinvgauss} uses a result from Chhikara and Folks (1974), with enhancements for right tails and log-probabilities.
\code{rinvgauss} uses an algorithm from Michael et al (1976).
\code{qinvgauss} uses code and algorithm from Giner and Smyth (2014).
}
 
\references{
Chhikara, R. S., and Folks, J. L., (1989).
\emph{The inverse Gaussian distribution: Theory, methodology and applications}.
Marcel Dekker, New York.

Chhikara, R. S., and Folks, J. L., (1974).
Estimation of the inverse Gaussian distribution function.
\emph{Journal of the American Statistical Association} 69, 250-254.

Giner, G., and Smyth, G. K. (2014).
A monotonically convergent Newton iteration for the quantiles of any unimodal distribution, with application to the inverse Gaussian distribution.
\url{http://www.statsci.org/smyth/pubs/qinvgaussPreprint.pdf}

Michael, J. R., Schucany, W. R., and Haas, R. W. (1976).
Generating random variates using transformations with multiple roots.
\emph{The American Statistician}, 30, 88--90.

Tweedie, M. C. (1957).
Statistical Properties of Inverse Gaussian Distributions I.
\emph{Annals of Mathematical Statistics} 28, 362-377.
}

\author{Gordon Smyth with contributions from Paul Bagshaw, Centre National d'Etudes des Telecommunications, France, and Trevor Park, Department of Statistics, University of Florida}

\seealso{
\code{dinvGauss}, \code{pinvGauss}, \code{qinvGauss} and \code{rinvGauss} in the SuppDists package.
}

\examples{
q <- rinvgauss(20,mean=1,dispersion=0.5) # generate vector of 20 random numbers
p <- pinvgauss(q,mean=1,dispersion=0.5) # p should be uniform
}

\keyword{distribution}
