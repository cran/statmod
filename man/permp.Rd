\name{permp}
\alias{permp}

\title{Exact permutation p-values}
\description{
Calculates exact p-values for permutation tests when permutations are randomly drawn with replacement.
}
\usage{
permp(x, nperm, n1, n2, total.nperm=NULL, method="auto", twosided=TRUE)
}
\arguments{
 \item{x}{number of permutations that yielded test statistics at least as extreme as the observed data. May be a vector or an array of values.}
 \item{nperm}{total number of permutations performed.}
 \item{n1}{sample size of group 1. Not required if \code{total.nperm} is supplied.}
 \item{n2}{sample size of group 2. Not required if \code{total.nperm} is supplied.}
 \item{total.nperm}{total number of permutations allowable from the design of the experiment.}
 \item{method}{character string indicating computation method. Possible values are \code{"exact"}, \code{"approximate"} or \code{"auto"}.}
 \item{twosided}{logical, is the test two-sided and symmetric between the two groups?}
}

\details{
This function can be used for calculating exact p-values for permutation tests where permutations are sampled with replacement, using theory and methods developed by Phipson and Smyth (2010).
The input values are the total number of permutations done (\code{nperm}) and the number of these that were considered at least as extreme as the observed data (\code{x}).

\code{total.nperm} is the total number of distinct values of the test statistic that are possible.
This is generally equal to the number of possible permutations, unless a two-sided test is conducted with equal sample sizes, in which case \code{total.nperm} is half the number of permutations, because the test statistic must then be symmetric in the two groups.
Usually \code{total.nperm} is computed automatically from \code{n1} and \code{n2}, but can also be supplied directly by the user.

When \code{method="exact"}, the p-values are computed to full machine precision by summing a series terms.
When \code{method="approximate"}, an approximation is used that is faster and uses less memory.
If \code{method="auto"}, the exact calculation is used when \code{total.nperm} is less than or equal to 10,000 and the approximation is used otherwise.
}
\value{
vector or array of p-values, of same dimensions as \code{x}
}
 
\author{Belinda Phipson and Gordon Smyth}

\references{
Phipson B, and Smyth GK (2010).
Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Issue 1, Article 39.
\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}
}

\examples{
x <- 0:5
# Both calls give same results
permp(x=x, nperm=99, n1=6, n2=6)
permp(x=x, nperm=99, total.nperm=462)
}

\keyword{htest}
