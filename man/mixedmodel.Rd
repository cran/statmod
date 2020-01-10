\name{mixedModel2}
\alias{mixedModel2}
\alias{mixedModel2Fit}
\alias{randomizedBlock}
\alias{randomizedBlockFit}

\title{Fit Mixed Linear Model with 2 Error Components}
\description{
Fits a mixed linear model by REML.  The linear model contains one random factor apart from the unit errors.
}

\usage{
mixedModel2(formula, random, weights=NULL, only.varcomp=FALSE, data=list(),
            subset=NULL, contrasts=NULL, tol=1e-6, maxit=50, trace=FALSE)
mixedModel2Fit(y, X, Z, w=NULL, only.varcomp=FALSE, tol=1e-6, maxit=50, trace=FALSE)
randomizedBlock(formula, random, weights=NULL, only.varcomp=FALSE, data=list(),
            subset=NULL, contrasts=NULL, tol=1e-6, maxit=50, trace=FALSE)
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
Note that \code{randomizedBlock} and \code{mixedModel2} are alternative names for the same function.

This function fits the model \eqn{y=Xb+Zu+e} where \eqn{b} is a vector of fixed coefficients and \eqn{u} is a vector of random effects.
Write \eqn{n} for the length of \eqn{y} and \eqn{q} for the length of \eqn{u}.
The random effect vector \eqn{u} is assumed to be normal, mean zero, with covariance matrix \eqn{\sigma^2_uI_q} while \eqn{e} is normal, mean zero, with covariance matrix \eqn{\sigma^2I_n}.
If \eqn{Z} is an indicator matrix, then this model corresponds to a randomized block experiment.
The model is fitted using an eigenvalue decomposition that transforms the problem into a Gamma generalized linear model.
To the knowledge of the author, this is an original and unpublished approach to the problem of fitting mixed models.

Note that the block variance component \code{varcomp[2]} is not constrained to be non-negative.
It may take negative values corresponding to negative intra-block correlations.
However the correlation \code{varcomp[2]/sum(varcomp)} must lie between \code{-1} and \code{1}.

Missing values in the data are not allowed.

This function is in principle equivalent to \code{lme(fixed=formula,random=~1|random)} except that the block variance component is not constrained to be non-negative.
If the block variance component is non-negative, then this function is faster and more accurate than \code{lme} for small to moderate size data sets but slower than \code{lme} when the number of observations is large.

This function tends to be fast and reliable, compared to competitor functions that fit randomized block models, when then number of observations is small, say no more than 200.
However it becomes quadratically slow as the number of observations increases because of the need to compute two singular value decompositions of order nearly equal to the number of observations, although this can be limited to only one decomposition if \code{only.varcomp = TRUE}).
For these reasons, this function is a good choice when fitting large numbers of mixed models to small datasets but is not optimal as currently implemented for fitting mixed models to very large data sets.
}

\value{
A list with the components:
\item{varcomp}{numeric vector of length two containing the residual and block components of variance.}
\item{se.varcomp}{standard errors for the variance components (if \code{only.varcomp=FALSE}).}
\item{coefficients}{numeric vector of fixed coefficients (if \code{only.varcomp=FALSE}).}
\item{se.coefficients}{standard errors for the fixed coefficients (if \code{only.varcomp=FALSE}).}
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
out <- mixedModel2(Y~SG+VP+V10+EP, random=No, data=petrol)
cbind(varcomp=out$varcomp,se=out$se.varcomp)
}

\keyword{regression}
