\name{glmnb.fit}
\alias{glmnb.fit}
\title{Fit Negative Binomial Generalized Linear Model with Log-Link}
\description{
Fit a generalized linear model with secure convergence.
}
\usage{
glmnb.fit(X, y, dispersion, weights = NULL, offset = 0, coef.start = NULL,
          start.method = "mean", tol = 1e-6, maxit = 50L, trace = FALSE)
}
\arguments{
  \item{X}{design matrix, assumed to be of full column rank.  Missing values not allowed.}
  \item{y}{numeric vector of responses. Negative or missing values not allowed.}
  \item{dispersion}{numeric vector of dispersion parameters for the negative binomial distribution. If of length 1, then the same dispersion is assumed for all observations.}
  \item{weights}{numeric vector of positive weights, defaults to all one.}
  \item{offset}{offset vector for linear model}
  \item{coef.start}{numeric vector of starting values for the regression coefficients}
  \item{start.method}{method used to find starting values, possible values are \code{"mean"} or \code{"log(y)"}}
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
This function implements a modified Fisher scoring algorithm for generalized linear models, analogous to the Levenberg-Marquardt algorithm for nonlinear least squares.
The Levenberg-Marquardt modification checks for a reduction in the deviance at each step, and avoids the possibility of divergence.
The result is a very secure algorithm that converges for almost all datasets.

\code{glmnb.fit} is in principle equivalent to \code{glm.fit(X,y,family=negative.binomial(link="log",theta=1/dispersion))} but with more secure convergence.
Here \code{negative.binomial} is a function in the MASS package.

The \code{dispersion} parameter is the same as \code{1/theta} for the \code{MASS::negative.binomial} function or \code{1/size} for the \code{stats::rnbinom} function.
\code{dispersion=0} corresponds to the Poisson distribution.
}
\author{Gordon Smyth and Yunshun Chen}
 
\references{
Dunn, PK, and Smyth, GK (2018).
\emph{Generalized linear models with examples in R}. Springer, New York, NY.
\doi{10.1007/978-1-4419-0118-7}
}

\seealso{
The \code{glmFit} function in the edgeR package on Bioconductor is a high-performance version of \code{glmnb.fit} for many \code{y} vectors at once.

\code{\link{glm}} is the standard glm fitting function in the stats package.
\code{negative.binomial} in the MASS package defines a negative binomial family for use with \code{glm}.
}
\examples{
y <- rnbinom(10, mu=1:10, size=5)
X <- cbind(1, 1:10)
fit <- glmnb.fit(X, y, dispersion=0.2, trace=TRUE)
}
\keyword{regression}
