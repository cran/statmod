\name{sage.test}
\alias{sage.test}

\title{Compare Two SAGE Libraries}
\description{
Compute a p-values for differential expression for each tag between two SAGE libraries. 
}

\usage{
sage.test(x, y, n1=sum(x), n2=sum(y))
}

\arguments{
\item{x}{integer vector giving counts in first library}
\item{y}{integer vector giving counts in second library}
\item{n1}{total number of tags in first library}
\item{n2}{total number of tags in second library}
}

\details{
This function uses a binomial approximation to the Fisher Exact test for each tag.
This will be accurate when the total number of tags is large.
}

\value{
Numeric vector of p-values.
}

\author{Gordon Smyth}

\seealso{
\code{\link{fisher.test2}}, \code{\link[ctest]{fisher.test}}
}

\examples{
sage.test(c(0,5,10),c(0,30,50),n1=10000,n2=15000)
}

\keyword{htest}
