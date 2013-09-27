\name{meanT}
\alias{meanT}

\title{Mean t-Statistic Between Two Groups of Growth Curves}
\description{
The mean-t statistic of the distance between two groups of growth curves. 
}

\usage{
meanT(y1, y2)
}

\arguments{
\item{y1}{matrix of response values for the first group, with a row for each individual and a column for each time. Missing values are allowed.}
\item{y2}{matrix of response values for the second group.  Must have the same number of columns as \code{y1}. Missing values are allowed.}
}

\details{
This function computes the pooled two-sample t-statistic between the response values at each time, and returns the mean of these values weighted by the degrees of freedom.
This function is used by \code{compareGrowthCurves}.
}

\value{numeric vector of length one containing the mean t-statistic.}

\author{Gordon Smyth}

\seealso{
\code{\link{compareGrowthCurves}}, \code{\link{compareTwoGrowthCurves}}
}

\examples{
y1 <- matrix(rnorm(4*3),4,3)
y2 <- matrix(rnorm(4*3),4,3)
meanT(y1,y2)

data(PlantGrowth)
compareGrowthCurves(PlantGrowth$group,as.matrix(PlantGrowth$weight))
# Can make p-values more accurate by nsim=10000
}

\keyword{regression}
