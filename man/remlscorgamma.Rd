\name{remlscoregamma}
\alias{remlscoregamma}
\title{Approximate REML for gamma regression with structured dispersion}
\description{
Estimates structured dispersion effects using approximate REML with gamma responses.
}
\usage{
remlscoregamma(y,X,Z,mlink="log",dlink="log",trace=FALSE,tol=1e-5,maxit=40)
}
\arguments{
  \item{y}{numeric vector of responses}
  \item{X}{design matrix for predicting the mean}
  \item{Z}{design matrix for predicting the variance}
  \item{mlink}{character string or numeric value specifying link for mean model}
  \item{dlink}{character string or numeric value specifying link for dispersion model}
  \item{trace}{Logical variable. If true then output diagnostic information at each iteration.}
  \item{tol}{Convergence tolerance}
  \item{maxit}{Maximum number of iterations allowed}
}
\value{
List with the following components:
  \item{beta}{Vector of regression coefficients for predicting the mean}
  \item{se.beta}{<Standard errors for beta}
  \item{gamma}{Vector of regression coefficients for predicting the variance}
  \item{se.gam}{Standard errors for gamma}
  \item{mu}{Estimated means}
  \item{phi}{Estimated dispersions}
  \item{deviance}{Minus twice the REML log-likelihood}
  \item{h}{Leverages}
}
\details{
Write \eqn{\mu_i=E(y_i)} for the expectation of the \eqn{i}th response and \eqn{s_i=\var(y_i)}.
We assume the heteroscedastic regression model
\deqn{\mu_i=\bold{x}_i^T\bold{\beta}}
\deqn{\log(\sigma^2_i)=\bold{z}_i^T\bold{\gamma},}
where \eqn{\bold{x}_i} and \eqn{\bold{z}_i} are vectors of covariates, and \eqn{\bold{\beta}} and \eqn{\bold{\gamma}} are vectors of regression coefficients affecting the mean and variance respectively.

Parameters are estimated by maximizing the REML likelihood using REML scoring as described in Smyth and Verbyla (2001).
See also Smyth and Verbyla (1999a,b).
}
\references{
Smyth, G. K., and Verbyla, A. P. (1999a). Adjusted likelihood methods for modelling dispersion in generalized linear models. \emph{Environmetrics} 10, 695-709. 
\url{http://www.statsci.org/smyth/pubs/ties98tr.html}

Smyth, G. K., and Verbyla, A. P. (1999b). Double generalized linear models: approximate REML and diagnostics. In \emph{Statistical Modelling: Proceedings of the 14th International Workshop on Statistical Modelling}, Graz, Austria, July 19-23, 1999, H. Friedl, A. Berghold, G. Kauermann (eds.), Technical University, Graz, Austria, pages 66-80.
\url{http://www.statsci.org/smyth/pubs/iwsm99-Preprint.pdf}

Smyth, GK, and Verbyla, AP (2009). Leverage adjustments for dispersion modelling in generalized nonlinear models. \emph{Australian and New Zealand Journal of Statistics} 51, 433-448.
}
\examples{
data(welding)
attach(welding)
y <- Strength
X <- cbind(1,(Drying+1)/2,(Material+1)/2)
colnames(X) <- c("1","B","C")
Z <- cbind(1,(Material+1)/2,(Method+1)/2,(Preheating+1)/2)
colnames(Z) <- c("1","C","H","I")
out <- remlscoregamma(y,X,Z)
}
\keyword{regression}
