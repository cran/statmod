\name{remlscoregamma}
\alias{remlscoregamma}
\title{Approximate REML for Gamma Regression with Structured Dispersion}
\description{
Estimates structured dispersion effects using approximate REML with gamma responses.
}
\usage{
remlscoregamma(y, X, Z, mlink="log", dlink="log", trace=FALSE, tol=1e-5, maxit=40)
}
\arguments{
  \item{y}{numeric vector of responses.}
  \item{X}{design matrix for predicting the mean.}
  \item{Z}{design matrix for predicting the variance.}
  \item{mlink}{character string or numeric value specifying link for mean model.}
  \item{dlink}{character string or numeric value specifying link for dispersion model.}
  \item{trace}{logical value. If \code{TRUE} then diagnostic information is output at each iteration.}
  \item{tol}{convergence tolerance.}
  \item{maxit}{maximum number of iterations allowed.}
}
\value{
List with the following components:
  \item{beta}{numeric vector of regression coefficients for predicting the mean.}
  \item{se.beta}{numeric vector of standard errors for beta.}
  \item{gamma}{numeric vector of regression coefficients for predicting the variance.}
  \item{se.gam}{numeric vector of standard errors for gamma.}
  \item{mu}{numeric vector of estimated means.}
  \item{phi}{numeric vector of estimated dispersions.}
  \item{deviance}{minus twice the REML log-likelihood.}
  \item{h}{numeric vector of leverages.}
}
\details{
This function fits a double generalized linear model (glm) with gamma responses.
As for ordinary gamma glms, a link-linear model is assumed for the expected values.
The double glm assumes a separate link-linear model for the dispersions as well.
The responses \code{y} are assumed to follow a gamma generalized linear model with link \code{mlink} and design matrix
\code{X}.
The dispersions follow a link-linear model with link \code{dlink} and design matrix \code{Z}.

Write \eqn{y_i} for the \eqn{i}th response.
The \eqn{y_i} are assumed to be independent and gamma distributed with \eqn{E(y_i) = \mu_i} and var\eqn{(y_i)=\phi_i\mu_i^2}.
The link-linear model for the means can be written as
\deqn{g(\mu)=X\beta}
where \eqn{g} is the mean-link function defined by \code{mlink} and \eqn{\mu} is the vector of means.
The dispersion link-linear model can be written as
\deqn{h(\phi)=Z\gamma}
where \eqn{h} is the dispersion-link function defined by \code{dlink} and \eqn{\phi} is the vector of dispersions.

The parameters \eqn{\gamma} are estimated by approximate REML likelihood using an adaption of the algorithm described by Smyth (2002).
See also Smyth and Verbyla (1999a,b) and Smyth and Verbyla (2009).
Having estimated \eqn{\gamma} and \eqn{\phi}, the \eqn{\beta} are estimated as usual for a gamma glm.

The estimated values for \eqn{\beta}, \eqn{\mu}, \eqn{\gamma} and \eqn{\phi} are return as \code{beta}, \code{mu}, \code{gamma} and \code{phi} respectively.
}
\references{
Smyth, G. K., and Verbyla, A. P. (1999a). Adjusted likelihood methods for modelling dispersion in generalized linear models. \emph{Environmetrics} 10, 695-709. 
\url{http://www.statsci.org/smyth/pubs/ties98tr.html}

Smyth, G. K., and Verbyla, A. P. (1999b). Double generalized linear models: approximate REML and diagnostics. In \emph{Statistical Modelling: Proceedings of the 14th International Workshop on Statistical Modelling}, Graz, Austria, July 19-23, 1999, H. Friedl, A. Berghold, G. Kauermann (eds.), Technical University, Graz, Austria, pages 66-80.
\url{http://www.statsci.org/smyth/pubs/iwsm99-Preprint.pdf}

Smyth, G. K. (2002). An efficient algorithm for REML in heteroscedastic regression. \emph{Journal of Computational and Graphical Statistics} \bold{11}, 836-847.

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
