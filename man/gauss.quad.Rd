\name{gauss.quad}
\alias{gauss.quad}
\title{Gaussian Quadrature}

\description{Calculate nodes and weights for Gaussian quadrature.}

\usage{gauss.quad(n,kind="legendre",alpha=0,beta=0)}

\arguments{
  \item{n}{number of nodes and weights}
  \item{kind}{kind of Gaussian quadrature, one of \code{"legendre"}, \code{"chebyshev1"}, \code{"chebyshev2"}, \code{"hermite"}, \code{"jacobi"} or \code{"laguerre"}}
  \item{alpha}{parameter for Jacobi or Laguerre quadrature, must be greater than -1}
  \item{beta}{parameter for Jacobi quadrature, must be greater than -1}
}

\value{A list containing the components
  \item{nodes}{vector of values at which to evaluate the function}
  \item{weights}{vector of weights to give the function values}
}

\details{
The integral from \code{a} to \code{b} of \code{w(x)*f(x)} is approximated by \code{sum(w*f(x))} where \code{x} is the vector of nodes and \code{w} is the vector of weights.  The approximation is exact if \code{f(x)} is a polynomial of order no more than \code{2n+1}.
The possible choices for \code{w(x)}, \code{a} and \code{b} are as follows:

Legendre quadrature: \code{w(x)=1} on \code{(-1,1)}.

Chebyshev quadrature of the 1st kind: \code{w(x)=1/sqrt(1-x^2)} on \code{(-1,1)}.

Chebyshev quadrature of the 2nd kind: \code{w(x)=sqrt(1-x^2)} on \code{(-1,1)}.

Hermite quadrature: \code{w(x)=exp(-x^2)} on \code{(-Inf,Inf)}.

Jacobi quadrature: \code{w(x)=(1-x)^alpha*(1+x)^beta} on \code{(-1,1)}. Note that Chebyshev quadrature is a special case of this.

Laguerre quadrature: \code{w(x)=x^alpha*exp(-x)} on \code{(0,Inf)}.

The method is explained in Golub and Welsch (1969).
}

\note{
This function solves a dense \code{n}x\code{n} eigenvector problem and is therefore slow for large \code{n}.
It could be made far more efficient by using an eigenvector function designed to compute the leading terms of the eigenvectors for tridiagonal matrices.
}

\references{
Golub, G. H., and Welsch, J. H. (1969). Calculation of Gaussian
quadrature rules. \emph{Mathematics of Computation} \bold{23}, 221-230.

Golub, G. H. (1973). Some modified matrix eigenvalue problems.
\emph{Siam Review} \bold{15}, 318-334.

Stroud and Secrest (1966). \emph{Gaussian Quadrature Formulas}. Prentice-
Hall, Englewood Cliffs, N.J.
}

\author{Gordon Smyth}

\seealso{
\code{\link[statmod:gauss.quad.prob]{gauss.quad.prob}}, \code{\link[base:integrate]{integrate}}
}

\examples{
out <- gauss.quad(10,"laguerre",alpha=5)
sum(out$weights * out$nodes) / gamma(6)
#  mean of gamma distribution with alpha=6
}
\keyword{math}
