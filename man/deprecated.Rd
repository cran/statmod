\name{deprecated}
\alias{statmod-deprecated}
\alias{randomizedBlock}
\alias{randomizedBlockFit}

\title{Deprecated Functions in statmod Package}
\description{
These functions are provided for compatibility with older versions
of R only, and may be defunct as soon as the next release.
}

\usage{
randomizedBlock(formula, random, weights=NULL, only.varcomp=FALSE, data=list(), subset=NULL, contrasts=NULL, tol=1e-6, maxit=50, trace=FALSE)
randomizedBlockFit(y, X, Z, w=NULL, only.varcomp=FALSE, tol=1e-6, maxit=50, trace=FALSE)
}

\arguments{
The arguments \code{formula}, \code{weights}, \code{data}, \code{subset} and \code{contrasts} have the same meaning as in \code{lm}.
The arguments \code{y}, \code{X} and \code{w} have the same meaning as in \code{lm.wfit}.
\item{formula}{formula specifying the fixed model.}
\item{random}{vector or factor specifying the blocks corresponding to random effects.}
\item{weights}{optional vector of prior weights.}
\item{only.varcomp}{logical value, if \code{TRUE} computation of standard errors and fixed effect coefficients will be skipped}
\item{data}{an optional data frame containing the variables in the model.}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{contrasts}{an optional list. See the \code{contrasts.arg} argument of \code{model.matrix.default}.}
\item{tol}{small positive numeric tolerance, passed to \code{glmgam.fit}}
\item{maxit}{maximum number of iterations permitted, passed to \code{glmgam.fit}}
\item{trace}{logical value, passed to \code{glmgam.fit}. If \code{TRUE} then working estimates will be printed at each iteration.}
\item{y}{numeric response vector}
\item{X}{numeric design matrix for fixed model}
\item{Z}{numeric design matrix for random effects}
\item{w}{optional vector of prior weights}
}

\details{
\code{randomizedBlock} and \code{randomizedBlockFit} have been renamed to \code{\link{mixedModel2}} and \code{\link{mixedModel2Fit}} respectively.
}

\keyword{regression}
