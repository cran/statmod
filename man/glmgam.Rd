\name{glmgam.fit}
\alias{glmgam.fit}
\title{Fit Gamma Generalized Linear Model by Fisher Scoring with Identity Link}
\description{
Fit a generalized linear model with secure convergence.
}
\usage{
glmgam.fit(X, y, coef.start = NULL, tol = 1e-6, maxit = 50, trace = FALSE)
}
\arguments{
  \item{X}{design matrix, assumed to be of full column rank.  Missing values not allowed.}
  \item{y}{numeric vector of responses. Negative or missing values not allowed.}
  \item{coef.start}{numeric vector of starting values for the regression coefficients}
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
This function implements a modified Fisher scoring algorithm for generalized linear models, similar to the Levenberg-Marquardt algorithm for nonlinear least squares.
The Levenberg-Marquardt modification checks for a reduction in the deviance at each step, and avoids the possibility of divergence.
The result is a very secure algorithm that converges for almost all datasets.

\code{glmgam.fit} is in principle equivalent to \code{glm.fit(X,y,family=Gamma(link="identity"))} but with much more secure convergence.
}

\author{Gordon Smyth and Yunshun Chen}
 
\references{
Dunn, PK, and Smyth, GK (2018).
\emph{Generalized linear models with examples in R}. Springer, New York, NY.
\doi{10.1007/978-1-4419-0118-7}
}

\seealso{
\code{glmgam.fit} is called by \code{\link{mixedModel2Fit}}.

\code{\link{glm}} is the standard glm fitting function in the stats package.
}

\examples{
y <- rgamma(10, shape=5)
X <- cbind(1, 1:10)
fit <- glmgam.fit(X, y, trace=TRUE)
}
\keyword{regression}
