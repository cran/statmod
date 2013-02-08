\name{plot.limdil}
\docType{class}
\alias{print.limdil}
\alias{plot.limdil}
\title{Plot or print an object of class limdil}

\description{
Plot or print the results of a limiting dilution analysis.
}

\usage{
\S3method{print}{limdil}(x, \dots)
\S3method{plot}{limdil}(x, col.group=NULL, cex=1, lwd=1, legend.pos="bottomleft", \dots)
}

\arguments{
\item{x}{object of class \code{limdil}.}
\item{col.group}{vector of colors for the groups of the same length as \code{levels(x$group)}.}
\item{cex}{relative symbol size}
\item{lwd}{relative line width}
\item{legend.pos}{positioning on plot of legend when there are multiple groups}
\item{\dots}{other arguments to be passed to \code{plot} or \code{print}.
Note that \code{pch} and \code{lty} are reserved arguments for the plot method.}
}


\details{
The print method formats results similarly to a regression or anova summary in R.

The plot method produces a plot of a limiting dilution experiment similar to that in Bonnefoix et al (2001). 
The basic design of the plot was made popular by Lefkovits and Waldmann (1979).

The plot shows frequencies and confidence intervals for the multiple groups. 
A novel feature is that assays with 100\% successes are included in the plot and are represented by down-pointing triangles. 
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
