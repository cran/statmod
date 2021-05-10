\name{qresiduals}
\alias{qresiduals}
\alias{qresid}
\alias{qres.binom}
\alias{qres.pois}
\alias{qres.nbinom}
\alias{qres.gamma}
\alias{qres.invgauss}
\alias{qres.tweedie}
\alias{qres.default}
\title{Randomized Quantile Residuals}

\description{
Compute randomized quantile residuals for generalized linear models.}

\usage{
qresiduals(glm.obj,dispersion=NULL)
qresid(glm.obj,dispersion=NULL)
qres.binom(glm.obj)
qres.pois(glm.obj)
qres.nbinom(glm.obj)
qres.gamma(glm.obj,dispersion=NULL)
qres.invgauss(glm.obj,dispersion=NULL)
qres.tweedie(glm.obj,dispersion=NULL)
qres.default(glm.obj,dispersion=NULL)
}

\arguments{
\item{glm.obj}{Object of class \code{glm}. The generalized linear model family is assumed to be binomial for \code{qres.binom}, poisson for
\code{qres.pois}, negative binomial for \code{qres.nbinom}, Gamma for
\code{qres.gamma}, inverse Gaussian for \code{qres.invgauss} or
tweedie for \code{qres.tweedie}.}

\item{dispersion}{a positive real number. Specifies the value of the
dispersion parameter for a Gamma or inverse Gaussian generalized linear
model if known. If \code{NULL}, the dispersion will be estimated by its
Pearson estimator.}
}

\value{Numeric vector of standard normal quantile residuals.}

\details{
Quantile residuals are based on the idea of inverting the
estimated distribution function for each observation to
obtain exactly standard normal residuals. In the case of
discrete distributions, such as the binomial and Poisson,
some randomization is introduced to produce continuous
normal residuals. Quantile residuals are the residuals of
choice for generalized linear models in large dispersion
situations when the deviance and Pearson residuals can be
grossly non-normal. Quantile residuals are the only
useful residuals for binomial or Poisson data when the
response takes on only a small number of distinct values.
}
 
\references{
Dunn, K. P., and Smyth, G. K. (1996). Randomized quantile residuals. \emph{Journal of Computational and Graphical Statistics} \bold{5}, 1-10.
\url{http://www.statsci.org/smyth/pubs/residual.html}

Dunn, PK, and Smyth, GK (2018).
\emph{Generalized linear models with examples in R}. Springer, New York, NY.
\doi{10.1007/978-1-4419-0118-7}
}

\author{Gordon Smyth}

\seealso{
\code{\link{residuals.glm}}
}

\examples{
#  Poisson example: quantile residuals show no granularity
y <- rpois(20,lambda=4)
x <- 1:20
fit <- glm(y~x, family=poisson)
qr <- qresiduals(fit)
qqnorm(qr)
abline(0,1)

#  Gamma example:
#  Quantile residuals are nearly normal while usual resids are not
y <- rchisq(20, df=1)
fit <- glm(y~1, family=Gamma)
qr <- qresiduals(fit, dispersion=2)
qqnorm(qr)
abline(0,1)

#  Negative binomial example:
if(require("MASS")) {
fit <- glm(Days~Age,family=negative.binomial(2),data=quine)
summary(qresiduals(fit))
fit <- glm.nb(Days~Age,link=log,data = quine)
summary(qresiduals(fit))
}
}

\keyword{regression}
