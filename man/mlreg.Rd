\name{mlreg}
\alias{mlreg}
\alias{mlreg.fit}
\alias{mlreg.fit.zero}

\title{Fit a Linear Model by Maximum Likelihood}
\description{
Fits a linear model by maximum likelihood assuming a variety of response
distributions. 
}

\usage{
mlreg.fit(X, y, weights=NULL, dist="logistic", init=NULL, scale=NULL)
mlreg.fit.zero(y, weights=NULL, dist="logistic", init=NULL, scale=NULL)
}

\arguments{
\item{X}{the design matrix. Rows containing missing values (in X or y) will be removed.}
\item{y}{numeric response vector. Missing values will be removed.}
\item{weights}{vector of non-negative weights.}
\item{dist}{character string giving the name of the response distribution. The possibilities are \code{"extreme"}, \code{"logistic"}, \code{"gaussian"}, \code{"weibull"}, \code{"exponential"}, \code{"rayleigh"},   
\code{"loggaussian"}, \code{"lognormal"}, \code{"loglogistic"} or \code{"t"} on 4 df. }
\item{init}{numeric vector of initial values for the parameters.}
\item{scale}{if specified then the scale parameter is fixed at the given value.}
}

\details{
This function is merely a convenient wrapper for calling the
\code{survreg.fit} function, which is part of the \code{survival} library by Terry Therneau.
It fits the model y = X*b + scale*e where b is the vector of regression coefficients and e is a vector of mean-zero errors, by maximum likelihood.

The function \code{mlreg.fit.zero} assumes that the mean is zero and fits y = scale*e, estimating only the scale parameter.
}

\value{See the documentation for \code{survreg.object}}

\author{Gordon Smyth}

\seealso{
\code{\link[survival:survreg]{survreg}}, \code{\link[survival:survreg.object]{survreg.object}}
}

\examples{
x <- 1:50
y <- x + 2*rnorm(50)
X <- cbind(1,x)
out <- mlreg.fit(X,y,dist="logistic")
}

\keyword{regression}
