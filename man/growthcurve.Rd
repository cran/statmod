\name{growthcurve}
\alias{compareGrowthCurves}
\alias{compareTwoGrowthCurves}
\alias{plotGrowthCurves}

\title{Compare Groups of Growth Curves}
\description{
Do all pairwise comparisons between groups of growth curves using a permutation test. 
}

\usage{
compareGrowthCurves(group, y, levels=NULL, nsim=100, fun=meanT, times=NULL,
                    verbose=TRUE, adjust="holm", n0=0.5)
compareTwoGrowthCurves(group, y, nsim=100, fun=meanT, n0=0.5)
plotGrowthCurves(group, y, levels=sort(unique(group)), times=NULL, col=NULL,...)
}

\arguments{
\item{group}{vector or factor indicating group membership.  Missing values are allowed in \code{compareGrowthCurves} but not in \code{compareTwoGrowthCurves}.}
\item{y}{matrix of response values with rows for individuals and columns for times.  The number of rows must agree with the length of \code{group}. Missing values are allowed.}
\item{levels}{a character vector containing the identifiers of the groups to be compared.  By default all groups with two more more members will be compared.}
\item{nsim}{number of permutations to estimated p-values.}
\item{fun}{a function defining the statistic used to measure the distance between two groups of growth curves.
Defaults to \code{\link{meanT}}.}
\item{times}{a numeric vector containing the column numbers on which the groups should be compared.
By default all the columns are used.}
\item{verbose}{should progress results be printed?}
\item{adjust}{method used to adjust for multiple testing, see \code{p.adjust}.}
\item{n0}{offset used for numerator and denominator of p-value calculation.}
\item{col}{vector of colors corresponding to distinct groups}
\item{...}{other arguments passed to \code{plot()}}
}

\details{
\code{compareTwoGrowthCurves} performs a permutation test of the difference between two groups of growth curves.
\code{compareGrowthCurves} does all pairwise comparisons between two or more groups of growth curves.

The permutation p-values are computed as p = (ngt + neq/2 + n0) / (nsim + n0) where ngt is the number of permutations with test statistics greater than observed, neq is the number of permuttation with test statistics equal to that observed, and n0 is an offset to avoid p-values of zero (Phipson & Smyth 2010).
The offset n0 improves the type I error rate control and can be interpreted as allowing for the observed data as one of the permutations.
High resolution p-values can be obtained by setting \code{nsim} to some large value, \code{nsim=10000} say.
}

\value{
\code{compareTwoGrowthCurves} returns a list with two components, \code{stat} and \code{p.value}, containing the observed statistics and the estimated p-value.  \code{compareGrowthCurves} returns a data frame with components
\item{Group1}{name of first group in a comparison}
\item{Group2}{name of second group in a comparison}
\item{Stat}{observed value of the statistic}
\item{P.Value}{permutation p-value}
\item{adj.P.Value}{p-value adjusted for multiple testing}
}

\author{Gordon Smyth}

\references{
Elso, C. M., Roberts, L. J., Smyth, G. K., Thomson, R. J., Baldwin, T. M., 
Foote, S. J., and Handman, E. (2004). Leishmaniasis host response loci 
(lmr13) modify disease severity through a Th1/Th2-independent pathway.  
\emph{Genes and Immunity} 5, 93-100.

Baldwin, T., Sakthianandeswaren, A., Curtis, J., Kumar, B., Smyth, G. K., Foote, S., and Handman, E. (2007).
Wound healing response is a major contributor to the severity of cutaneous leishmaniasis in the ear model of infection.
\emph{Parasite Immunology} 29, 501-513.

Phipson B, Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Issue 1, Article 39.
\doi{10.2202/1544-6115.1585},
\doi{10.48550/arXiv.1603.05766}.
}

\seealso{
\code{\link{meanT}}, \code{\link{compareGrowthCurves}}, \code{\link{compareTwoGrowthCurves}}
}

\examples{
# A example with only one time
data(PlantGrowth)
compareGrowthCurves(PlantGrowth$group,as.matrix(PlantGrowth$weight))
# Can make p-values more accurate by nsim=10000
}

\keyword{regression}
