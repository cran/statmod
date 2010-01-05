\name{elda}
\alias{elda}
\alias{limdil}
\alias{limdil.class}

\title{Extreme Limiting Dilution Analysis}
\description{
Fit single-hit model to a dilution series using complementary log-log binomial regression.
}

\usage{
elda(response, dose, tested=rep(1,length(response)), group=rep(1,length(response)), observed=FALSE, confidence=0.95, test.unit.slope=FALSE)
limdil(response, dose, tested=rep(1,length(response)), group=rep(1,length(response)), observed=FALSE, confidence=0.95, test.unit.slope=FALSE)
}

\arguments{
\item{response}{numeric of integer counts of positive cases, out of \code{tested} trials}
\item{dose}{numeric vector of expected number of cells in assay}
\item{tested}{numeric vector giving number of trials at each dose}
\item{group}{vector or factor giving group to which the response belongs}
\item{observed}{logical, is the actual number of cells observed?}
\item{confidence}{numeric level for confidence interval}
\item{test.unit.slope}{logical, should the adequacy of the single-hit model be tested?}
}

\details{
This function is an implementation of maximum likelihood analysis of limiting dilution data with added features to accommodate small sample sizes (Hu and Smyth, 2009).
In particular, the function accommodates gracefully situations where 0% or 100% of the assays give positive results.
The methodology has typically been applied to the analysis of stem cell assays (Shackleton et al, 2006).

\code{elda} and \code{limdil} are alternative names for the same function.

A binomial generalized linear model is fitted for each group with cloglog link and offset \code{log(dose)}.
If \code{observed=FALSE}, a classic Poisson single-hit model is assumed, and the Poisson frequency of the stem cells is the \code{exp} of the intercept.
If \code{observed=TRUE}, the values of \code{dose} are treated as actual cell numbers rather than expected values.
This doesn't changed the generalized linear model fit but changes how the frequencies are extracted from the estimated model coefficient.

The confidence interval is a Wald confidence interval, unless all the responses are zero or at the maximum value, in which case Clopper-Pearson intervals are computed.

If \code{group} takes several values, then separate confidence intervals are computed for each group.
In this case it also possible to test for non-equality in frequency between the groups.
}

\value{
\code{limdil} produces an object of class \code{limdil} with the following components. There are \code{\link[=print.limdil]{print}} and \code{\link[=plot.limdil]{plot}} methods for limdil objects.

\item{CI}{numeric vector giving estimated frequency and lower and upper limits of Wald confidence interval of each group}
\item{test.difference}{numeric vector giving chisquare likelihood ratio test statistic and p-value for testing the difference between groups}
\item{test.slope.wald}{numeric vector giving wald test statistics and p-value for testing the slope of the offset equal to one}
\item{test.slope.lr}{numeric vector giving chisquare likelihood ratio test statistics and p-value for testing the slope of the offset equal to one}
\item{test.slope.scorel}{numeric vector giving score test statistics and p-value for testing multi-hit alternatives}
\item{test.slope.score}{numeric vector giving score test statistics and p-value for testing heterogeneity}
\item{response}{numeric of integer counts of positive cases, out of \code{tested} trials}
\item{tested}{numeric vector giving number of trials at each dose}
\item{dose}{numeric vector of expected number of cells in assay}
\item{group}{vector or factor giving group to which the response belongs}
\item{num.group}{number of groups}
}

\author{Yifang Hu and Gordon Smyth}

\references{
Shackleton, M., Vaillant, F., Simpson, K. J., Stingl, J., Smyth, G. K., Asselin-Labat, M.-L., Wu, L., Lindeman, G. J., and Visvader, J. E. (2006).
Generation of a functional mammary gland from a single stem cell.
\emph{Nature} 439, 84-88.
\url{http://www.nature.com/nature/journal/v439/n7072/abs/nature04372.html}

Hu, Y, and Smyth, GK (2009).
ELDA: Extreme limiting dilution analysis for comparing depleted and enriched populations in stem cell and other assays.
\emph{Journal of Immunological Methods} 347, 70-78.
\url{http://dx.doi.org/10.1016/j.jim.2009.06.008}
}

\examples{
# When there is one group
Dose <- c(50,100,200,400,800)
Responses <- c(2,6,9,15,21)
Tested <- c(24,24,24,24,24)
out <- limdil(Responses,Dose,Tested,test.unit.slope=TRUE)
out
plot(out)

# When there are four groups
Dose <- c(30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500)
Responses <- c(2,3,2,1,6,5,6,1,2,3,4,2,6,6,6,1)
Tested <- c(6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6)
Group <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)
}

\keyword{regression}
