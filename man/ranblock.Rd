\name{ranblock}
\alias{randomizedBlock}
\alias{randomizedBlockFit}

\title{Randomized Block Mixed Linear Model}
\description{
Fits a mixed linear model by REML.  The linear model must contain only one random factor apart from the unit errors.
}

\usage{
randomizedBlock(formula, random, weights=NULL, fixed.estimates=TRUE, data=list(), subset=NULL, contrasts=NULL)
randomizedBlockFit(y,X,Z,w=NULL,fixed.estimates=TRUE)
}

\arguments{
The arguments \code{formula}, \code{weights}, \code{data}, \code{subset} and \code{contrasts} have the same meaning as in \code{lm}.
The arguments \code{y}, \code{X} and \code{w} have the same meaning as in \code{lm.wfit}.
\item{formula}{formula specifying the fixed model.}
\item{random}{vector or factor specifying the blocks corresponding to random effects.}
\item{weights}{optional vector of prior weights.}
\item{fixed.estimates}{should the fixed effect coefficients be returned?}
\item{data}{an optional data frame containing the variables in the model.}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.}
\item{y}{response vector.}
\item{X}{design matrix for fixed model.}
\item{Z}{design matrix for random effects.}
\item{w}{optional vector of prior weights.}
}

\details{
This function fits the model \eqn{y=Xb+Zu+e} where \eqn{b} is a vector of fixed coefficients and \eqn{u} is a vector of random effects.
Write \eqn{n} for the length of \eqn{y} and \eqn{q} for the length of \eqn{u}.
The random effect vector \eqn{u} is assumed to be normal, mean zero, with covariance matrix \eqn{\sigma^2_uI_q} while \eqn{e} is normal, mean zero, with covariance matrix \eqn{\sigma^2I_n}.
If \eqn{Z} is an indicator matrix, then this model corresponds to a randomized block experiment.
The model is fitted using an eigenvalue decomposition which transforms the problem into a Gamma generalized linear model.

This function is essentially equivalent to \code{lme(fixed=formula,random=~1|random)} but is more accurate and is much faster for small to moderate size data sets.

Note that the block variance component \code{sigmasquared[2]} is not constrained to be non-negative.
It may take negative values corresponding to negative intra-block correlations.
However the correlation \code{sigmasquared[2]/sum(sigmasquared)} should lie between \code{-1} and \code{1}.

Missing values in the data are not allowed.
}

\value{A list with the components.
If \code{fixed.estimates=TRUE} then the components from \code{"lm.fit"} are also returned.
\item{sigmasquared}{vector of length two containing the residual and block components of variance.}
\item{se.sigmasquared}{standard errors for the components of variance.}
}

\author{Gordon Smyth}

\seealso{
\code{\link{glmgam.fit}}, \code{\link[nlme:lme]{lme}}, \code{\link[base:lm]{lm}}, \code{\link[base:lm]{lm.fit}}
}

\examples{
#  Compare with first data example from Venable and Ripley (2002), Chapter 10, "Linear Models"
library(MASS)
data(petrol)
out <- randomizedBlock(Y~SG+VP+V10+EP, random=No, data=petrol)
cbind(sigmasquared=out$sigmasquared,se=out$se.sigmasquared)
}

\keyword{regression}
