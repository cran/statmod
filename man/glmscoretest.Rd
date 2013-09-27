\name{glm.scoretest}
\alias{glm.scoretest}

\title{Score Test for Adding a Covariate to a GLM}
\description{
Computes score test statistics (z-statistics) for adding covariates to a generalized linear model.
}

\usage{
glm.scoretest(fit, x2, dispersion=NULL)
}

\arguments{
\item{fit}{generalized linear model fit object, of class \code{\link{glm}}.}
\item{x2}{vector or matrix with each column a covariate to be added.}
\item{dispersion}{the dispersion for the generalized linear model family.}
}

\details{
Rao's score statistic.
Is the locally most powerful test for testing vs a one-sided alternative.
Asympotically equivalent to likelihood ratio tests, but convenient for one-sided tests.

This function computes a score test statistics for adding each covariate individually.

The dispersion parameter is treated as for \code{\link{summary.glm}}.
If \code{NULL}, the Pearson estimator is used, except for the binomial and Poisson
families, for which the dispersion is one.
}

\value{numeric vector containing the z-statistics, one for each covariate.}

\author{Gordon Smyth}

\seealso{
\code{\link{glm}}, \code{\link{add1}}
}

\references{
Lovison, G (2005). On Rao score and Pearson $X^2$ statistics in generalized linear models.
\emph{Statistical Papers}, 46, 555-574.

Pregibon, D (1982). Score tests in GLIM with applications.
In \emph{GLIM82: Proceedings of the International Conference on Generalized Linear Models},
R Gilchrist (ed.), Lecture Notes in Statistics, Volume 14, Springer, New York, pages 87-97. 

Smyth, G. K. (2003). Pearson's goodness of fit statistic as a score test statistic. In: \emph{Science and Statistics: A Festschrift for Terry Speed}, D. R. Goldstein (ed.), IMS Lecture Notes - Monograph Series, Volume 40, Institute of Mathematical Statistics, Beachwood, Ohio, pages 115-126.
\url{http://www.statsci.org/smyth/pubs/goodness.pdf}
}

\examples{
#  Pearson's chisquare test for independence
#  in a contingency table is a score test.

#  First the usual test

y <- c(20,40,40,30)
chisq.test(matrix(y,2,2),correct=FALSE)

#  Now same test using glm.scoretest

a <- gl(2,1,4)
b <- gl(2,2,4)
fit <- glm(y~a+b,family=poisson)
x2 <- c(0,0,0,1)
z <- glm.scoretest(fit,x2)
z^2
}

\keyword{regression}
