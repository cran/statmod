\name{limdil}
\alias{limdil}

\title{Limiting Dilution Analysis}
\description{
Fit single-hit model to a dilution series using complementary log-log binomial regression.
}

\usage{
limdil(response,dose,tested=rep(1,length(response)),observed=FALSE,confidence=0.95,test.unit.slope=FALSE)
}

\arguments{
\item{response}{numeric of integer counts of positive cases, out of \code{tested} trials}
\item{dose}{numeric vector of expected number of cells in assay}
\item{tested}{numeric vector giving number of trials at each dose}
\item{observed}{logical, is the actual number of cells observed?}
\item{confidence}{numeric level for confidence interval}
\item{test.unit.slope}{logical, should the adequacy of the single-hit model be tested?}
}

\details{
A binomial generalized linear model is fitted with cloglog link and offset \code{log(dose)}.
If \code{observed=FALSE}, a classic Poisson single-hit model is assumed, and the Poisson frequency of the stem cells is the \code{exp} of the intercept.
If \code{observed=TRUE}, the values of \code{dose} are treated as actual cell numbers rather than expected values.
This doesn't changed the generalized linear model fit but changes how the frequencies are extracted from the estimated model coefficient.
}

\value{
List with components
\item{CI}{numeric vector giving estimated frequency and lower and upper limits of Wald confidence interval}
\item{test.unit.slope}{numeric vector giving chisquare likelihood ratio test statistic and p-value for testing the slope of the offset equal to one}
}

\author{Gordon Smyth}

\references{
Bonnefoix T, Bonnefoix P, Verdiel P, Sotto JJ. (1996).
Fitting limiting dilution experiments with generalized linear models results in a test of the single-hit Poisson assumption.
\emph{J Immunol Methods} 194, 113-119.
}

\examples{
Dose <- c(50,100,200,400,800)
Responses <- c(2,6,9,15,21)
Tested <- c(24,24,24,24,24)
limdil(Responses,Dose,Tested,test.unit.slope=TRUE)
}

\keyword{regression}
