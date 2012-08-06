\name{glmgam.fit}
\alias{glmgam.fit}
\alias{glmnb.fit}
\title{Fit Generalized Linear Model by Fisher Scoring with Levenberg Damping}
\description{
Fit a generalized linear model with secure convergence.
Provided for gamma glm with identity links or negative binomial glm with log-links.
}
\usage{
glmgam.fit(X, y, coef.start=NULL, tol=1e-6, maxit=50, trace=FALSE)
glmnb.fit(X, y, dispersion, offset=0, coef.start=NULL, start.method="mean", tol=1e-6, maxit=50, trace=FALSE)
}
\arguments{
  \item{X}{design matrix, assumed to be of full column rank.  Missing values not allowed.}
  \item{y}{numeric vector of responses. Negative or missing values not allowed.}
  \item{dispersion}{numeric vector of over-dispersion parameters for negative binomial. If of length 1, then same over-dispersion is assumed for all observations.}
  \item{offset}{offset vector for linear model}
  \item{coef.start}{numeric vector of starting values for the regression coefficients}
  \item{start.method}{method used to find starting values, possible values are \code{"mean"} and \code{"log(y)"}}
  \item{tol}{small positive numeric value giving convergence tolerance}
  \item{maxit}{maximum number of iterations allowed}
  \item{trace}{logical value. If \code{TRUE} then output diagnostic information at each iteration.}
}
\value{
List with the following components:
  \item{coefficients}{numeric vector of regression coefficients}
  \item{fitted}{numeric vector of fitted values}
  \item{deviance}{residual deviance}
  \item{iter}{number of iterations used to convergence. If convergence was not achieved then \code{iter} is set to \code{maxit+1}.}
}
\details{
These functions implement a modified Fisher scoring algorithm for generalized linear models, similar to the Levenberg-Marquardt algorithm for nonlinear least squares.
The Levenberg-Marquardt modification checks for a reduction in the deviance at each step, and avoids the possibility of divergence.
The result is a very secure algorithm that converges for almost all datasets.

\code{glmgam.fit} is in principle similar to \code{glm.fit(X,y,family=Gamma(link="identity"))} but with much more secure convergence.
This function is used by \code{\link{mixedModel2Fit}}.

\code{glmnb.fit} is in principle similar to \code{glm.fit(X,y,family=negative.binomial(link="log",theta=1/dispersion))} but with more secure convergence.

}
\author{Gordon Smyth and Yunshun Chen}
\examples{
y <- rgamma(10,shape=5)
X <- cbind(1,1:10)
fit <- glmgam.fit(X,y,trace=TRUE)
}
\keyword{regression}
