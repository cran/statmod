\name{growthcurve}
\alias{compareGrowthCurves}
\alias{compareTwoGrowthCurves}

\title{Compare Groups of Growth Curves}
\description{
Do all pairwise comparisons between groups of growth curves using a permutation test. 
}

\usage{
compareGrowthCurves(group,y,levels=NULL,nsim=100,fun=meanT,times=NULL,verbose=TRUE,adjust="holm")
compareTwoGrowthCurves(group,y,nsim=100,fun=meanT)
}

\arguments{
\item{group}{vector or factor indicating group membership.  Missing values are allowed in \code{compareGrowthCurves} but not in \code{compareTwoGrowthCurves}.}
\item{y}{matrix of response values with rows for individuals and columns for times.  The number of rows must agree with the length of \code{group}. Missing values are allowed.}
\item{levels}{a character vector containing the identifiers of the groups to be compared.  By default all groups with two more more members will be compared.}
\item{nsim}{number of permutations to estimated p-values.}
\item{fun}{the statistic used to measure the distance between two groups of growth curves.
Default to a mean t-statistic.}
\item{times}{a numeric vector containing the column numbers on which the groups should be compared.
By default all the columns are used.}
\item{verbose}{should progress results be printed?}
\item{adjust}{method used to adjust for multiple testing, see \code{p.adjust}.}
}

\details{
\code{compareTwoGrowthCurves} performs a permutation test of the difference between two groups of growth curves.
\code{compareGrowthCurves} does all pairwise comparisons between two or more groups of growth curves.
Accurate p-values can be obtained by setting \code{nsim} to some large value, \code{nsim=10000} say.
}

\value{
\code{compareTwoGrowthCurves} returns a list with two components, \code{stat} and \code{p.value}, containing the observed statistics and the estimated p-value.  \code{compareGrowthCurves} returns a data frame with components
\item{Group1}{name of first group in a comparison}
\item{Group2}{name of second group in a comparison}
\item{Stat}{observed value of the statistic}
\item{P.Value}{estimated p-value}
\item{adj.P.Value}{p-value adjusted for multiple testing}
}

\author{Gordon Smyth}

\seealso{
\code{\link{compareGrowthCurves}}, \code{\link{compareTwoGrowthCurves}}
}

\examples{
# A example with only one time
data(PlantGrowth)
compareGrowthCurves(PlantGrowth$group,as.matrix(PlantGrowth$weight))
# Can make p-values more accurate by nsim=10000
}

\keyword{regression}
