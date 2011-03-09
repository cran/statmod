\name{sage.test}
\alias{sage.test}

\title{Exact Binomial Tests For Comparing Two SAGE Libraries}
\description{
Compute binomial test p-values for differential expression for each tag between two serial analysis of gene expression (SAGE) libraries. 
}

\usage{
sage.test(x, y, n1=sum(x), n2=sum(y))
}

\arguments{
\item{x}{integer vector giving counts in first library.
Non-integer values are rounded to the nearest integer.}
\item{y}{integer vector giving counts in second library.
Non-integer values are rounded to the nearest integer.}
\item{n1}{total number of tags in first library.
Non-integer values are rounded to the nearest integer.}
\item{n2}{total number of tags in second library.
Non-integer values are rounded to the nearest integer.}
}

\details{
SAGE is a method for counting the frequency of sequence tags in samples of RNA.
One can test for tagwise differential expression between the SAGE libraries for two different RNA samples using an exact binomial test, conditioning on the total count for each tag.
When the total counts is very large, a Pearson chisquare approximation is used.

The binomial test is similar to Fisher's exact test, but probably even more appropriate in this context.
}

\value{
Numeric vector of p-values.
}

\author{Gordon Smyth}

\references{
\url{http://www.sagenet.org}
}

\seealso{
\code{\link{binom.test}}
}

\examples{
sage.test(c(0,5,10),c(0,30,50),n1=10000,n2=15000)
#  Univariate equivalents:
binom.test(5,5+30,p=10000/(10000+15000))$p.value
binom.test(10,10+50,p=10000/(10000+15000))$p.value
}

\keyword{htest}
