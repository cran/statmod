\name{sage.test}
\alias{sage.test}

\title{Compare Two SAGE Libraries}
\description{
Compute p-values for differential expression for each tag between two SAGE libraries. 
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
This function uses a binomial approximation to the Fisher Exact test for each tag.
The approximation is accurate when \code{n1} and \code{n2} are large and \code{x} and \code{y} are small in comparison.
}

\value{
Numeric vector of p-values.
}

\author{Gordon Smyth}

\seealso{
\code{\link[ctest]{fisher.test}}
}

\examples{
sage.test(c(0,5,10),c(0,30,50),n1=10000,n2=15000)
#  Exact equivalents
fisher.test(matrix(c(0,0,10000-0,15000-0),2,2))$p.value
fisher.test(matrix(c(5,30,10000-5,15000-30),2,2))$p.value
fisher.test(matrix(c(10,50,10000-10,15000-50),2,2))$p.value
}

\keyword{htest}
