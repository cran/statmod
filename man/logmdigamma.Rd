\name{logmdigamma}
\alias{logmdigamma}
\title{Log Minus Digamma Function}
\description{
The difference between the \code{log} and \code{digamma} functions.
}
\usage{
logmdigamma(x)
}
\arguments{
  \item{x}{numeric vector or array of positive values. Negative or zero values will return \code{NA}.}
}
\details{
\code{digamma(x)} is asymptotically equivalent to \code{log(x)}.  \code{logmdigamma(x)} computes \code{log(x) - digamma(x)} without subtractive cancellation for large \code{x}.
}
\author{Gordon Smyth}
\references{
Abramowitz, M., and Stegun, I. A. (1970). \emph{Handbook of
mathematical functions.} Dover, New York.
}
\seealso{
\code{\link{digamma}}
}
\examples{
log(10^15) - digamma(10^15) # returns 0
logmdigamma(10^15) # returns value correct to 15 figures
}
\keyword{math}
