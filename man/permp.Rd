\name{permp}
\alias{permp}

\title{Exact permutation p-values}
\description{
Calculates exact p-values for permutation tests with permutations sampled with replacement.
}
\usage{
permp(x, nperm, n1, n2, total.nperm=NULL, method="auto", twosided=TRUE)
}
\arguments{
 \item{x}{numeric vector or array giving the number of permutations that yield test statistics at least as extreme as that observed.}
 \item{nperm}{number of permutations performed.}
 \item{n1}{sample size of group 1. Not required if \code{total.nperm} is supplied.}
 \item{n2}{sample size of group 2. Not required if \code{total.nperm} is supplied.}
 \item{total.nperm}{total number of permutations allowable from the design of the experiment.}
 \item{method}{computation method, possible values are \code{"exact"}, \code{"approximate"} or \code{"auto"}.}
 \item{twosided}{logical, is the test two-sided so that the test statistic is symmetric between the two groups?}
}

\details{
This function can be used for calculating exact p-values for permutation tests where permutations are sampled with replacement, using theory and methods developed by Phipson and Smyth (2010).
\code{total.nperm} is the total number of distinct values of the test statistic that are possible.
This is generally equal to the number of possible permutations, unless a two-sided test is conducted with equal sample sizes, in which case \code{total.nperm} is half the number of permutations, because the test statistic must then be symmetric in the two groups.
\code{total.nperm} can be supplied directly by the user, alternatively it will be computed from \code{n1} and \code{n2}.

When \code{method="exact"}, the p-values are computed to full machine precision by evaluating a summation.
When \code{method="approximate"}, an approximation is used that is faster and uses less memory.
If \code{method="auto"}, the exact calculation is used when \code{total.nperm} is less than or equal to 10,000, otherwise the approximation is used.
}
\value{
vector or array of p-values, of same dimensions as \code{x}
}
 
\author{Belinda Phipson and Gordon Smyth}

\references{
Phipson, B and Smyth, GK (2010). Exact permutation p-values when permutations are randomly drawn. Technical Report, Walter and Eliza Hall Institute of Medical Research, Melbourne, Australia.
}

\examples{
#  Consider a two-sided permutation test with 99 permutations, 
#  5 reaching the statistically significant threshold. 
#  Assume a two group experiment with 6 in each group.
#  Input total.nperm=462

permp(x=5, nperm=99, total.nperm=462)

# Input n1=6 and n2=6

permp(x=5, nperm=99, n1=6, n2=6)

# Suppose we have a vector of p-values taking on the values 0 to 10 for the same experiment described above.

x<-0:10
permp(x=x, nperm=99, total.nperm=462)
}

\keyword{htest}
