\name{mscale}
\alias{mscale}

\title{M Scale Estimation}
\description{
Robust estimation of a scale parameter using Hampel's redescending psi function.
}

\usage{
mscale(u, na.rm=FALSE)
}

\arguments{
\item{u}{numeric vector of residuals.}
\item{na.rm}{logical. Should missing values be removed?}
}

\value{numeric constant giving the estimated scale.}

\details{
Estimates a scale parameter or standard deviation using an M-estimator with 50\% breakdown.
This means the estimator is highly robust to outliers.
If the input residuals \code{u} are a normal sample, then \code{mscale(u)} should be equal to the standard deviation.
}

\author{Gordon Smyth}

\references{
Yohai, V. J. (1987). High breakdown point and high efficiency robust estimates for regression. \emph{Ann. Statist.} 15, 642-656.

Stromberg, A. J. (1993). Computation of high breakdown nonlinear regression parameters. \emph{J. Amer. Statist. Assoc.} 88, 237-244.

Smyth, G. K., and Hawkins, D. M. (2000). Robust frequency estimation using elemental sets. \emph{Journal of Computational and Graphical Statistics} 9, 196-214. 
}

%\seealso{
%\code{\link{rho.hampel}}, \code{\link{psi.hampel}}
%}

\examples{
u <- rnorm(100)
sd(u)
mscale(u)
}
