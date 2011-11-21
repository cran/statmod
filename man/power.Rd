\name{power.fisher.test}
\alias{power.fisher.test}
\title{Power of Fisher's Exact Test for Comparing Proportions}
\description{
Calculate by simulation the power of Fisher's exact test for comparing two proportions given two margin counts. 
}

\usage{
power.fisher.test(p1, p2, n1, n2, alpha=0.05, nsim=100, alternative="two.sided")
}

\arguments{
\item{p1}{first proportion to be compared.}
\item{p2}{second proportion to be compared.}
\item{n1}{first sample size.}
\item{n2}{second sample size.}
\item{alpha}{significance level.}
\item{nsim}{number of data sets to simulate.}
\item{alternative}{indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".}
}

\details{
Computes the power of Fisher's exact test for testing the null hypothesis that \code{p1} equals \code{p2} against the alternative that they are not equal.
}

\value{
Estimated power of the test.
}

\author{Gordon Smyth}

\seealso{
\code{\link{fisher.test}}, \code{\link{power.t.test}}
}

\examples{
power.fisher.test(0.5,0.9,20,20) # 70% chance of detecting difference
}

\keyword{htest}
