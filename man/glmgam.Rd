\name{glmgam.fit}
\alias{glmgam.fit}
\title{Gamma Generalized Linear Model with Identity Link}
\description{
Estimates a gamma generalized linear model with identity link using Fisher scoring with Levenberg damping.
}
\usage{
glmgam.fit(X,y,start=NULL,trace=FALSE,tol=1e-6,maxit=50)
}
\arguments{
  \item{X}{design matrix, assumed to be of full column rank.  Missing values not allowed.}
  \item{y}{numeric vector of responses. Missing values not allowed.}
  \item{start}{numeric vector of starting values for the regression coefficients}
  \item{trace}{Logical variable. If true then output diagnostic information at each iteration.}
  \item{tol}{Convergence tolerance}
  \item{maxit}{Maximum number of iterations allowed}
}
\value{
List with the following components:
  \item{coefficients}{numeric vector of regression coefficients}
  \item{fitted}{numeric vector of fitted values}
  \item{deviance}{residual deviance}
}
\item{details}
This function is similar to \code{glm.fit(X,y,family=Gamma(link="identity"))} but has more secure convergence.

This function is used by \code{randomizedBlockFit}.
}
\examples{
y <- rgamma(10,shape=5)
X <- cbind(1,1:10)
fit <- glmgam.fit(X,y,trace=TRUE)
}
\keyword{regression}
