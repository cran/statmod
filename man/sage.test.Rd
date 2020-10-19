\name{sage.test}
\alias{sage.test}

\title{Exact Binomial Tests For Comparing Two SAGE Libraries (Obsolete)}
\description{
Computes p-values for differential abundance for each tag between two digital libraries,
conditioning on the total count for each tag.
The counts in each group as a proportion of the whole are assumed to follow a binomial distribution.
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
This function was originally written for comparing SAGE libraries (a method for counting the frequency of sequence tags in samples of RNA).
It can however be used for comparing any two digital libraries from RNA-Seq, ChIP-Seq or other technologies with respect to technical variation.

An exact two-sided binomial test is computed for each tag.
This test is closely related to Fisher's exact test for 2x2 contingency tables but, unlike Fisher's test, it conditions on the total number of counts for each tag.
The null hypothesis is that the expected counts are in the same proportions as the library sizes, i.e., that the binomial probability for the first library is \code{n1/(n1+n2)}.

The two-sided rejection region is chosen analogously to Fisher's test.
Specifically, the rejection region consists of those values with smallest probabilities
under the null hypothesis.

When the counts are reasonably large, the binomial test, Fisher's test and Pearson's chisquare all give the same results.
When the counts are smaller, the binomial test is usually to be preferred in this context.

This function is a later version of the earlier \code{sage.test} function in the sagenhaft Bioconductor package.
This function has been made obsolete by \code{binomTest} in the edgeR package.
}

\note{
This function is kept in the statmod package so as not to break code that depends on it
but it has been replaced by \code{binomTest} in the edgeR Bioconductor package and 
is no longer updated.
It may be removed in a later release of this package.
}

\value{
Numeric vector of p-values.
}

\references{
\url{https://en.wikipedia.org/wiki/Binomial_test}

\url{https://en.wikipedia.org/wiki/Fisher's_exact_test}

\url{https://en.wikipedia.org/wiki/Serial_analysis_of_gene_expression}

\url{https://en.wikipedia.org/wiki/RNA-Seq}
}

\author{Gordon Smyth}

\seealso{
The \code{binomTest} function in the \href{https://doi.org/doi:10.18129/B9.bioc.edgeR}{edgeR} package on Bioconductor is a newer and better version of this function.

\code{\link{binom.test}} in the stats package performs univariate binomial tests.
}

\examples{
sage.test(c(0,5,10),c(0,30,50),n1=10000,n2=15000)
#  Univariate equivalents:
binom.test(5,5+30,p=10000/(10000+15000))$p.value
binom.test(10,10+50,p=10000/(10000+15000))$p.value
}

\keyword{htest}
