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
Write \eqn{\mu_i=E(y_i)} for the expectation of the $i$th response and $s_i=\var(y_i)$.
We assume the heteroscedastic regression model
\deqn{\mu_i=\bold{x}_i^T\bold{\beta}}
\deqn{\log(\sigma^2_i)=\bold{z}_i^T\bold{\gamma},}
where $\bold{x}_i$ and $\bold{z}_i$ are vectors of covariates, and $\bold{\beta}$ and $\bold{\gamma}$ are vectors of regression coefficients affecting the mean and variance respectively.

Parameters are estimated by maximizing the REML likelihood using REML scoring as described in Smyth (2002).
}
\references{
Smyth, G. K., and Verbyla, A. P. (2002). Leverage adjustments for dispersion modelling in generalized nonlinear models.
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
