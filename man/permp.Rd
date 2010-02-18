\name{permp}
\alias{permp}

\title{Exact permutation p-values}
\description{
Calculates exact p-values for permutation tests with permutations sampled with replacement.
}
\usage{
permp(x, nperm, n1, n2, total.nperm=NULL)
}
\arguments{
 \item{x}{the number of cases that reach the significance threshold in the permutation test}
 \item{nperm}{the number of permutations performed}
 \item{n1}{sample size of group 1}
 \item{n2}{sample size of group 2}
 \item{total.nperm}{the total number of permutations allowable from the design of the experiment}
  }
\details{
This function can be used for calculating p-values for permutation tests where permutations are sampled with replacement. If the total number of permutations, total.nperm, is not known, it can be calculated by specifying the sample sizes in each group, n1 and n2. Note that this is only for a two group experiment. This function calculates the exact permutation p-value, which is given by: (x+1)/(m+1) - (Integration Term). The integration term uses total.nperm, and is enumerated using Gaussian quadrature. It works for a scalar or vector of x's.
}
\value{
	\code{permp} outputs either a single p-value or a vector of p-values depending on the input
 }
 
\author{Belinda Phipson and Gordon Smyth}
\keyword{htest}

\examples{
#  Consider a permutation test with 99 permutations, 
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

