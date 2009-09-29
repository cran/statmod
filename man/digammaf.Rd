\name{Digamma}
\alias{Digamma}
\alias{canonic.digamma}
\alias{d2cumulant.digamma}
\alias{unitdeviance.digamma}
\alias{cumulant.digamma}
\alias{meanval.digamma}
\alias{varfun.digamma}
\title{Digamma generalized linear model family}
\description{
Produces a Digamma generalized linear model family object. The Digamma distribution is the
distribution of the unit deviance for a gamma response.
}
\usage{
Digamma(link = "log")
unitdeviance.digamma(y, mu)
cumulant.digamma(theta)
meanval.digamma(theta)
d2cumulant.digamma(theta)
varfun.digamma(mu)
canonic.digamma(mu)
}
\arguments{
  \item{link}{character string, number or expressing specifying the link function. See \code{quasi} for specification of this argument.}
  \item{y}{numeric vector of (positive) response values}
  \item{mu}{numeric vector of (positive) fitted values}
  \item{theta}{numeric vector of values of the canonical variable, equal to \eqn{-1/\phi} where \eqn{\phi} is the dispersion parameter of the gamma distribution}
}
\value{
\code{Digamma} produces a glm family object, which is a list of functions and expressions used by \code{glm} in its iteratively reweighted least-squares algorithm. See \code{family} for details.

The other functions take vector arguments and produce vector values of the same length and called by \code{Digamma}.
\code{unitdeviance.digamma} gives the unit deviances of the family, equal to the squared deviance residuals.
\code{cumulant.digamma} is the cumulant function.  If the dispersion is unity, then successive derivatives of the cumulant function give successive cumulants of the Digamma distribution.  \code{meanvalue.digamma} gives the first derivative, which is the expected value.
\code{d2cumulant.digamma} gives the second derivative, which is the variance.
\code{canonic.digamma} is the inverse of \code{meanvalue.digamma} and gives the canonical parameter as a function of the mean parameter.
\code{varfun.digamma} is the variance function of the Digamma family, the variance as a function of the mean.
}
\details{
This family is useful for dispersion modelling with gamma generalized linear models.
The Digamma distribution describes the distribution of the unit deviances for a gamma family, in the same way that the gamma distribution itself describes the distribution of the unit deviances for Gaussian or inverse Gaussian families.
The Digamma distribution is so named because it is dual to the gamma distribution in the above sense, and because the \code{digamma function} appears in its mean function.

Suppose that \eqn{y} follows a gamma distribution with mean \eqn{\mu} and dispersion parameter \eqn{\phi}, so the variance of \eqn{y} is \eqn{\phi \mu^2}.
Write \eqn{d(y,\mu)} for the gamma distribution unit deviance.
Then \code{meanval.digamma(-1/phi)} gives the mean of \eqn{d(y,\mu)} and \code{2*d2cumulant.digamma(-1/phi)} gives the variance.
}

\author{Gordon Smyth}

\references{
Smyth, G. K. (1989). Generalized linear models with varying dispersion. \emph{J. R. Statist. Soc. B}, \bold{51}, 47-61.
}
\examples{
# Test for log-linear dispersion trend in gamma regression
y <- rchisq(20,df=1)
x <- 1:20
out.gam <- glm(y~x,family=Gamma(link="log"))
d <- residuals(out.gam)^2
out.dig <- glm(d~x,family=Digamma(link="log"))
summary(out.dig,dispersion=2)
}
\seealso{
\code{\link{quasi}}, \code{\link{make.link}}
}
\keyword{models}
