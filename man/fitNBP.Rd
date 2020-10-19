\name{fitNBP}
\alias{fitNBP}

\title{Negative Binomial Model for SAGE Libraries with Pearson Estimation of Dispersion}
\description{
Fit a multi-group negative-binomial model to SAGE data, with Pearson estimation of the common overdispersion parameter. 
}

\usage{
fitNBP(y, group=NULL, lib.size=colSums(y), tol=1e-5, maxit=40, verbose=FALSE)
}

\arguments{
\item{y}{numeric matrix giving counts. Rows correspond to tags (genes) and columns to SAGE libraries.}
\item{group}{factor indicating which library belongs to each group. If \code{NULL} then one group is assumed.}
\item{lib.size}{vector giving total number of tags in each library.}
\item{tol}{small positive numeric tolerance to judge convergence}
\item{maxit}{maximum number of iterations permitted}
\item{verbose}{logical, if \code{TRUE} then iteration progress information is output.}
}

\details{
The overdispersion parameter is estimated equating the Pearson goodness of fit to its expectation.
The variance is assumed to be of the form
Var(y)=mu*(1+phi*mu)
where E(y)=mu and phi is the dispersion parameter.
All tags are assumed to share the same dispersion.

For given dispersion, the model for each tag is a negative-binomial generalized linear model with log-link and \code{log(lib.size)} as offset.
The coefficient parametrization used is that corresponding to the formula \code{~0+group+offset(log(lib.size)}.

Except for the dispersion being common rather than genewise, the model fitted by this function is equivalent to that proposed by Lu et al (2005).
The numeric algorithm used is that of alternating iterations (Smyth, 1996) using Newton's method as the outer iteration for the dispersion parameter starting at phi=0.
This iteration is monotonically convergent for the dispersion.
}

\note{
This function has been made obsolete by the \href{https://doi.org/doi:10.18129/B9.bioc.edgeR}{edgeR} package on Bioconductor.
}

\value{
List with components
\item{coefficients}{numeric matrix of rates for each tag (gene) and each group}
\item{fitted.values}{numeric matrix of fitted values}
\item{dispersion}{estimated dispersion parameter}
}

\author{Gordon Smyth}

\references{
Lu, J, Tomfohr, JK, Kepler, TB (2005).
Identifying differential expression in multiple SAGE libraries: an overdispersed log-linear model approach.
\emph{BMC Bioinformatics} 6,165.

Smyth, G. K. (1996). Partitioned algorithms for maximum likelihood and other nonlinear estimation. \emph{Statistics and Computing}, 6, 201-216.
}

\seealso{
\code{\link{sage.test}}

The edgeR package on Biconductor provides new and better functions to fit negative-binomial glms to SAGE or RNA-seq data.
See particularly the \code{glmFit} and \code{mglmOneWay} functions.
}

\examples{
# True value for dispersion is 1/size=2/3
# Note the Pearson method tends to under-estimate the dispersion
y <- matrix(rnbinom(10*4,mu=4,size=1.5),10,4)
lib.size <- rep(50000,4)
group <- c(1,1,2,2)
fit <- fitNBP(y,group=group,lib.size=lib.size)
logratio <- fit$coef \%*\% c(-1,1)
}

\keyword{regression}
