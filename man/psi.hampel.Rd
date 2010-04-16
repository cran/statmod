\name{psi.hampel}
\alias{psi.hampel}

\title{Hampel's Redescending psi Function}
\description{
Hampel's redescending psi function.
}

\usage{
psi.hampel(u, a = 1.5, b = 3.5, c = 8)
}

\arguments{
\item{u}{numeric vector of residuals.}
\item{a}{numeric constant.}
\item{b}{numeric constant.}
\item{c}{numeric constant.}
}

\value{numeric vector like \code{u} giving the values of Hampel's psi function.}

\details{
psi is equal to \code{u} up to \code{a}, then is flat to \code{b}, then decreases linearly to zero at \code{c}.
}

\author{Gordon Smyth}

\references{
Hampel, Ronchetti, Rousseeuw and Stahel (1986). \emph{Robust Statistics.} Wiley, New York, page 150.
	
Yohai, V. J. (1987). High breakdown point and high efficiency robust estimates for regression. \emph{Ann. Statist.} 15, 642-656.
	
Stromberg, A. J. (1993). Computation of high breakdown nonlinear regression parameters. \emph{J. Amer. Statist. Assoc.} 88, 237-244.
	
Smyth, G. K., and Hawkins, D. M. (2000). Robust frequency estimation using elemental sets. \emph{Journal of Computational and Graphical Statistics} 9, 196-214.
}

\seealso{
\code{\link{mscale}}, \code{\link{rho.hampel}}
}