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
dinvgauss(x, mu, lambda=1)
pinvgauss(q, mu, lambda=1)
qinvgauss(p, mu, lambda=1)
rinvgauss(n, mu, lambda=1)
}

\arguments{
  \item{x}{vector of quantiles. Missing values (NAs) are allowed.}
  \item{q}{vector of quantiles. Missing values (NAs) are allowed.}
  \item{p}{vector of probabilities. Missing values (NAs) are allowed.}
  \item{n}{sample size. If \code{length(n)} is larger than 1, then \code{length(n)} random values are returned.}
  \item{mu}{vector of (positive) means. This is replicated to be the same length as \code{p} or \code{q} or the number of deviates generated.}
  \item{lambda}{vector of (positive) precision parameters. This is replicated to be the same length as \code{p} or \code{q} or the number of deviates generated.}
}

\value{
Vector of same length as \code{x} or \code{q} giving the density (\code{dinvgauss}), probability (\code{pinvgauss}), quantile (\code{qinvgauss}) or random sample (\code{rinvgauss}) for the inverse
Gaussian distribution with mean \code{mu} and inverse dispersion \code{lambda}.
Elements of \code{q} or \code{p} that are missing will cause the corresponding elements of
the result to be missing.
}

\details{
The inverse Gaussian distribution takes values on the positive real line. The variance of the distribution is $\mu^3/\lambda$. Applications of the inverse Gaussian include sequential analysis, diffusion processes and radiotechniques. The inverse Gaussian is one of the response distributions used in generalized linear models.
}
 
\references{
Chhikara, R. S., and Folks, J. Leroy, (1989). \emph{The inverse Gaussian distribution: Theory, methodology, and  applications}. Marcel Dekker, New York.
}

\author{Gordon Smyth and Dr Paul Bagshaw, Centre National d'Etudes des Telecommunications (DIH/DIPS), France (\code{qinvgauss})}

\seealso{
\code{\link{qres}}, \code{\link[base]{inverse.gaussian}}, \code{\link[base]{.Random.seed}} about random number generation.
}

\examples{
y <- rinvgauss(20,1,2) # generate vector of 20 random numbers
p <- pinvgauss(y,1,2) # p should be uniform</tt>
}

\keyword{distribution}
