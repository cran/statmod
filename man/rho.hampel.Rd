\name{rho.hampel}
\alias{rho.hampel}

\title{Hampel's Redescending rho Function}
\description{
Integral of Hampel's redescending psi function.
}

\usage{
rho.hampel(u, a = 1.5, b = 3.5, c = 8)
}

\arguments{
\item{u}{numeric vector of residuals.}
\item{a}{numeric constant.}
\item{b}{numeric constant.}
\item{c}{numeric constant.}
}

\value{numeric vector like u giving the values of Hampel's rho function.}

\details{
The values of rho are analogous sum of squared residuals. rho is actually the integral of \code{\link[=psi.hampel]{Hampel's redescending psi function}}. The default values for the parameters are as in Stromberg (1993).
}

\author{Gordon Smyth}

\references{
Hampel, Ronchetti, Rousseeuw and Stahel (1986). \emph{Robust Statistics.} Wiley, New York, page 150.

Yohai, V. J. (1987). High breakdown point and high efficiency robust estimates for regression. \emph{Ann. Statist.} 15, 642-656.

Stromberg, A. J. (1993). Computation of high breakdown nonlinear regression parameters. \emph{J. Amer. Statist. Assoc.} 88, 237-244.

Smyth, G. K., and Hawkins, D. M. (2000). Robust frequency estimation using elemental sets. \emph{Journal of Computational and Graphical Statistics} 9, 196-214.
}

\seealso{
\code{\link{mscale}}, \code{\link{psi.hampel}}
}
