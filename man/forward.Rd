\name{forward}
\alias{forward}

\title{Forward Selection of Covariates for Multiple Regression}
\description{
Fit a multi-group negative-binomial model to SAGE data, with Pearson estimation of the common overdispersion parameter. 
}

\usage{
forward(y, x, xkept=NULL, intercept=TRUE, nvar=ncol(x))
}

\arguments{
\item{y}{numeric response vector.}
\item{x}{numeric matrix of covariates, candidates to be added to the regression.}
\item{xkept}{numeric matrix of covariates to be included in the starting regression.}
\item{intercept}{logical, should an intercept be added to \code{xkept}?}
\item{nvar}{integer, number of covariates from \code{x} to add to the regression.}
}

\details{
This function has the advantage that \code{x} can have many more columns than the length of \code{y}.
}

\value{
Integer vector of length \code{nvar}, giving the order in which columns of \code{x} are added to the regression.
}

\author{Gordon Smyth}

\seealso{
\code{\link{step}}
}

\examples{
y <- rnorm(10)
x <- matrix(rnorm(10*5),10,5)
forward(y,x)
}

\keyword{regression}
