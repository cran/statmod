\name{elda}
\alias{elda}
\alias{limdil}
\alias{eldaOneGroup}
\alias{limdil.class}
\alias{limdil-class}

\title{Extreme Limiting Dilution Analysis}
\description{
Fit single-hit model to a dilution series using complementary log-log binomial regression.
}

\usage{
elda(response, dose, tested=rep(1,length(response)), group=rep(1,length(response)),
     observed=FALSE, confidence=0.95, test.unit.slope=FALSE)
limdil(response, dose, tested=rep(1,length(response)), group=rep(1,length(response)),
     observed=FALSE, confidence=0.95, test.unit.slope=FALSE)
eldaOneGroup(response, dose, tested, observed=FALSE, confidence=0.95,
     tol=1e-8, maxit=100, trace=FALSE)
}

\arguments{
\item{response}{numeric vector giving number of positive cases out of \code{tested} trials. Should take non-negative integer values.}
\item{dose}{numeric vector of expected number of cells in assay. Values must be positive.}
\item{tested}{numeric vector giving number of trials at each dose. Should take integer values.}
\item{group}{vector or factor giving group to which the response belongs.}
\item{observed}{logical, is the actual number of cells observed?}
\item{confidence}{numeric level for confidence interval. Should be strictly between 0 and 1.}
\item{test.unit.slope}{logical, should the adequacy of the single-hit model be tested?}
\item{tol}{convergence tolerance.}
\item{maxit}{maximum number of Newton iterations to perform.}
\item{trace}{logical, if \code{TRUE} then iterim results are output at each iteration.}
}

\details{
\code{elda} and \code{limdil} are alternative names for the same function.
(\code{limdil} was the older name before the 2009 paper by Hu and Smyth.)
\code{eldaOneGroup} is a lower-level function that does the computations when there is just one group, using a globally convergent Newton iteration.
It is called by the other functions.

These functions implement maximum likelihood analysis of limiting dilution data using methods proposed by Hu and Smyth (2009).
The functions gracefully accommodate situations where 0\% or 100\% of the assays give positive results, which is why we call it "extreme" limiting dilution analysis.
The functions provide the ability to test for differences in stem cell frequencies between groups, and to test goodness of fit in a number of ways.
The methodology has been applied to the analysis of stem cell assays (Shackleton et al, 2006).

The statistical method is explained by Hu and Smyth (2009).
A binomial generalized linear model is fitted for each group with cloglog link and offset \code{log(dose)}.
If \code{observed=FALSE}, a classic Poisson single-hit model is assumed, and the Poisson frequency of the stem cells is the \code{exp} of the intercept.
If \code{observed=TRUE}, the values of \code{dose} are treated as actual cell numbers rather than expected values.
This doesn't change the generalized linear model fit, but it does change how the frequencies are extracted from the estimated model coefficient (Hu and Smyth, 2009).

The confidence interval is a Wald confidence interval, unless the responses are all negative or all positive, in which case Clopper-Pearson intervals are computed.

If \code{group} takes several values, then separate confidence intervals are computed for each group.
In this case a likelihood ratio test is conducted for differences in active cell frequencies between the groups.

These functions compute a number of different tests of goodness of fit.
One test is based on the coefficient for \code{log(dose)} in the generalized linear model.
The nominal slope is 1.
A slope greater than one suggests a multi-hit model in which two or more cells are synergistically required to produce a positive response.
A slope less than 1 suggests some sort of cell interference.
Slopes less than 1 can also be due to heterogeneity of the stem cell frequency between assays.
\code{elda} conducts likelihood ratio and score tests that the slope is equal to one.

Another test is based on the coefficient for \code{dose}.
This idea is motivated by a suggestion of Gart and Weiss (1967), who suggest that heterogeneity effects are more likely to be linear in \code{dose} than \code{log(dose)}.
These functions conducts score tests that the coefficient for \code{dose} is non-zero.
Negative values for this test suggest heterogeneity.

These functions produce objects of class \code{"limdil"}.
There are \code{\link[=print.limdil]{print}} and \code{\link[=plot.limdil]{plot}} methods for \code{"limdil"} objects.
}

\value{
\code{elda} and \code{limdil} produce an object of class \code{"limdil"}.  This is a list with the following components:
\item{CI}{numeric matrix giving estimated stem cell frequency and lower and upper limits of Wald confidence interval for each group}
\item{test.difference}{numeric vector giving chisquare likelihood ratio test statistic and p-value for testing the difference between groups}
\item{test.slope.wald}{numeric vector giving wald test statistics and p-value for testing the slope of the offset equal to one}
\item{test.slope.lr}{numeric vector giving chisquare likelihood ratio test statistics and p-value for testing the slope of the offset equal to one}
\item{test.slope.score.logdose}{numeric vector giving score test statistics and p-value for testing multi-hit alternatives}
\item{test.slope.score.dose}{numeric vector giving score test statistics and p-value for testing heterogeneity}
\item{response}{numeric of integer counts of positive cases, out of \code{tested} trials}
\item{tested}{numeric vector giving number of trials at each dose}
\item{dose}{numeric vector of expected number of cells in assay}
\item{group}{vector or factor giving group to which the response belongs}
\item{num.group}{number of groups}
}

\author{Yifang Hu and Gordon Smyth}

\references{
Hu, Y, and Smyth, GK (2009).
ELDA: Extreme limiting dilution analysis for comparing depleted and enriched populations in stem cell and other assays.
\emph{Journal of Immunological Methods} 347, 70-78.
\doi{10.1016/j.jim.2009.06.008}
\url{http://www.statsci.org/smyth/pubs/ELDAPreprint.pdf}

Shackleton, M., Vaillant, F., Simpson, K. J., Stingl, J., Smyth, G. K., Asselin-Labat, M.-L., Wu, L., Lindeman, G. J., and Visvader, J. E. (2006).
Generation of a functional mammary gland from a single stem cell.
\emph{Nature} 439, 84-88.
\url{https://www.nature.com/articles/nature04372}

Gart, JJ, and Weiss, GH (1967).
Graphically oriented tests for host variability in dilution experiments.
\emph{Biometrics} 23, 269-284.
}

\seealso{
\code{\link{plot.limdil}} and \code{\link{print.limdil}} are methods for \code{limdil} class objects.

A web interface to this function is available at \url{http://bioinf.wehi.edu.au/software/elda/}.
}

\examples{
# When there is one group
Dose <- c(50,100,200,400,800)
Responses <- c(2,6,9,15,21)
Tested <- c(24,24,24,24,24)
out <- elda(Responses,Dose,Tested,test.unit.slope=TRUE)
out
plot(out)

# When there are four groups
Dose <- c(30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500)
Responses <- c(2,3,2,1,6,5,6,1,2,3,4,2,6,6,6,1)
Tested <- c(6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6)
Group <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
elda(Responses,Dose,Tested,Group,test.unit.slope=TRUE)
}

\keyword{regression}
