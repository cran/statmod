\name{fisher.test2}
\alias{fisher.test2}

\title{Fast Version of Fisher's Exact Test for Count Data}

\description{
Performs Fisher's exact test for testing the null of independence of rows and columns in a contingency table with fixed marginals.
This is a cut down version of \code{fisher.test} in the \code{ctest} package.
}

\usage{
fisher.test2(x, workspace=2e+05)
}

\arguments{
\item{x}{numeric matrix containing contingency table counts}
\item{workspace}{an integer specifying the size of the workspace used in the network algorithm}
}

\value{
Numeric vector of p-values.
}

\author{Gordon Smyth}

\seealso{
\code{\link[ctest]{fisher.test}}
}

\examples{
fisher.test2(matrix(c(0,100,20,100),2,2))
}

\keyword{htest}
