\name{ranblock}
\alias{randomizedBlock}
\alias{randomizedBlockFit}

\title{Randomized Block Mixed Linear Model}
\description{
Fits a mixed linear model by REML.  The linear model must contain only one random factor apart from the unit errors.
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
This function fits the model \eqn{y=Xb+Zu+e} where \eqn{b} is a vector of fixed coefficients and \eqn{u} is a vector of random effects.
Write \eqn{n} for the length of \eqn{y} and \eqn{q} for the length of \eqn{u}.
The random effect vector \eqn{u} is assumed to be normal, mean zero, with covariance matrix \eqn{\sigma^2_uI_q} while \eqn{e} is normal, mean zero, with covariance matrix \eqn{\sigma^2I_n}.
If \eqn{Z} is an indicator matrix, then this model corresponds to a randomized block experiment.
The model is fitted using an eigenvalue decomposition which transforms the problem into a Gamma generalized linear model.

Note that the block variance component \code{varcomp[2]} is not constrained to be non-negative.
It may take negative values corresponding to negative intra-block correlations.
However the correlation \code{varcomp[2]/sum(varcomp)} must lie between \code{-1} and \code{1}.

Missing values in the data are not allowed.

This function is equivalent to \code{lme(fixed=formula,random=~1|random)}, except that the block variance component is not constrained to be non-negative, but is faster and more accurate for small to moderate size data sets.
It is slower than \code{lme} when the number of observations is large.

This function tends to be fast and reliable, compared to competitor functions which fit randomized block models, when then number of observations is small, say no more than 200.
However it becomes quadratically slow as the number of observations increases because of the need to do two eigenvalue decompositions of order nearly equal to the number of observations.
So it is a good choice when fitting large numbers of small data sets, but not a good choice for fitting large data sets.
}

\value{
A list with the components:
\item{varcomp}{vector of length two containing the residual and block components of variance.}
\item{se.varcomp}{standard errors for the components of variance.}
If \code{fixed.estimates=TRUE} then the components from \code{lm.fit} are also returned.
}

\author{Gordon Smyth}

\references{
Venables, W., and Ripley, B. (2002). \emph{Modern Applied Statistics with S-Plus}, Springer.
}

\seealso{
\code{\link{glmgam.fit}}, \code{\link[nlme:lme]{lme}}, \code{\link{lm}}, \code{\link{lm.fit}}
}

\examples{
#  Compare with first data example from Venable and Ripley (2002),
#  Chapter 10, "Linear Models"
library(MASS)
data(petrol)
out <- randomizedBlock(Y~SG+VP+V10+EP, random=No, data=petrol)
cbind(varcomp=out$varcomp,se=out$se.varcomp)
}

\keyword{regression}
