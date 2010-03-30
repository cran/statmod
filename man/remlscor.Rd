\name{remlscore}
\alias{remlscore}
\title{REML for Heteroscedastic Regression}
\description{
Fits a heteroscedastic regression model using residual maximum likelihood (REML).
}
\usage{
remlscore(y, X, Z, trace=FALSE, tol=1e-5, maxit=40)
}
\arguments{
  \item{y}{numeric vector of responses}
  \item{X}{design matrix for predicting the mean}
  \item{Z}{design matrix for predicting the variance}
  \item{trace}{Logical variable. If true then output diagnostic information at each iteration.}
  \item{tol}{Convergence tolerance}
  \item{maxit}{Maximum number of iterations allowed}
}
\value{
List with the following components:
  \item{beta}{vector of regression coefficients for predicting the mean}
  \item{se.beta}{vector of standard errors for beta}
  \item{gamma}{vector of regression coefficients for predicting the variance}
  \item{se.gam}{vector of standard errors for gamma}
  \item{mu}{estimated means}
  \item{phi}{estimated variances}
  \item{deviance}{minus twice the REML log-likelihood}
  \item{h}{numeric vector of leverages}
  \item{cov.beta}{estimated covariance matrix for beta}
  \item{cov.gam}{estimated covarate matrix for gamma}
  \item{iter}{number of iterations used}
}
\details{
Write \eqn{\mu_i=E(y_i)} for the expectation of the \eqn{i}th response and \eqn{s_i=\var(y_i)}.
We assume the heteroscedastic regression model
\deqn{\mu_i=\bold{x}_i^T\bold{\beta}}
\deqn{\log(\sigma^2_i)=\bold{z}_i^T\bold{\gamma},}
where \eqn{\bold{x}_i} and \eqn{\bold{z}_i} are vectors of covariates, and \eqn{\bold{\beta}} and \eqn{\bold{\gamma}} are vectors of regression coefficients affecting the mean and variance respectively.

Parameters are estimated by maximizing the REML likelihood using REML scoring as described in Smyth (2002).
}
\references{
Smyth, G. K. (2002). An efficient algorithm for REML in heteroscedastic regression. \emph{Journal of Computational and Graphical Statistics} \bold{11}, 836-847.
}
\author{Gordon Smyth}
\examples{
data(welding)
attach(welding)
y <- Strength
# Reproduce results from Table 1 of Smyth (2002)
X <- cbind(1,(Drying+1)/2,(Material+1)/2)
colnames(X) <- c("1","B","C")
Z <- cbind(1,(Material+1)/2,(Method+1)/2,(Preheating+1)/2)
colnames(Z) <- c("1","C","H","I")
out <- remlscore(y,X,Z)
cbind(Estimate=out$gamma,SE=out$se.gam)
}
\keyword{regression}
