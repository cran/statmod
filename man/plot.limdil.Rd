\name{plot.limdil}
\docType{class}
\alias{print.limdil}
\alias{plot.limdil}
\title{Plot or print an object of class limdil}

\description{
Plot or print an object of class limdil.
}

\usage{
\S3method{print}{limdil}(x, \dots)
\S3method{plot}{limdil}(x, col.group=NULL, \dots)
}

\arguments{
\item{x}{object of class \code{limdil}.}
\item{\dots}{other arguments to be passed to \code{plot} or \code{print}.
Note that \code{pch} and \code{lty} are reserved arguments for the plot method.}
\item{col.group}{vector of colors for the groups of the same length as \code{levels(x$group)}.}
}


\details{
Produces a plot of a limiting dilution experiment similar to that in Bonnefoix et al. (2001). 
The basic design of the plot was made popular by Lefkovits and Waldmann (1979).

The plot shows frequencies and confidence intervals for the multiple groups. 
A novel feature is that assays with 100\% successes are represented on the plot, by down-pointing triangles. 
}

\author{Yifang Hu and Gordon Smyth}

\references{
Bonnefoix, T, Bonnefoix, P, Callanan, M, Verdiel, P, and Sotto, JJ (2001).
Graphical representation of a generalized linear model-based statistical test estimating the fit of the single-hit poisson model to limiting
dilution assays.
\emph{The Journal of Immunology} 167, 5725-5730.

Lefkovits, I, and Waldmann, H (1979).
\emph{Limiting dilution analysis of cells in the immune system}.
Cambridge University Press, Cambridge.
}

\seealso{\link{limdil} describes the \code{limdil} class.}

\keyword{regression}