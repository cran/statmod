\name{hommel.test}
\alias{hommel.test}
\title{Test Multiple Comparisons}

\description{Given a set of p-values and a test level, returns test results for each of the hypotheses.}

\usage{hommel.test(p, alpha=0.05)
}

\arguments{
  \item{p}{vector of p-values}
  \item{alpha}{desired significance level}
}

\value{A logical vector indicating whether each hypothesis is accepted.}

\details{
This function implements the multiple testing procedure of Hommel (1988).
Hommel's method is also implemented as an adjusted P-value method in the function \code{p.adjust} but the accept/reject approach in this function is faster.
}

\references{
Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. \emph{Biometrika}, \bold{75}, 383-386.

Shaffer, J. P. (1995). Multiple hypothesis testing. \emph{Annual Review of Psychology} \bold{46}, 561-576. (An excellent review of the area.)
}

\author{Gordon Smyth}

\seealso{
\code{\link[statmod:p.adjust]{p.adjust}}
}

\examples{
p <- sort(runif(100))[1:10]
cbind(p,p.adjust(p,"hommel"),hommel.test(p))
}
\keyword{htest}
