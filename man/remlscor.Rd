\name{remlscore}
\alias{remlscore}
\title{REML for heteroscedastic regression}
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
  \item{beta}{Vector of regression coefficients for predicting the mean}
  \item{se.beta}{<Standard errors for beta}
  \item{gamma}{Vector of regression coefficients for predicting the variance}
  \item{se.gam}{Standard errors for gamma}
  \item{mu}{Estimated means}
  \item{phi}{Estimated variances}
  \item{deviance}{Minus twice the REML log-likelihood}
  \item{h}{Leverages}
}
\item{details}
Write \eqn{\mu_i=E(y_i)} for the expectation of the $i$th response and $s_i=\var(y_i)$.
We assume the heteroscedastic regression model
\deqn{\mm_i=\bold{x}_i^T\bold{\beta}}
\deqn{\log(\sigma^2_i=\bold{z}_i^T\bold{\gamma},}
where $\bold{x}_i$ and $\bold{z}_i$ are vectors of covariates, and $\bold{\beta}$ and $\bold{\gamma}$ are vectors of regression coefficients affecting the mean and variance respectively.

Parameters are estimated by maximizing the REML likelihood using REML scoring as described in Smyth (2002).
}
\references{
Smyth, G. K. (2002). An efficient algorithm for REML in heteroscedastic regression. \emph{Journal of Computational and Graphical Statistics} \bold{11}, 1-12.
}
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
