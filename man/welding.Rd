\name{welding}
\alias{welding}
\title{Data: Tensile Strength of Welds}
\description{
This is a highly fractionated two-level factorial design employed as a screening design
in an off-line welding experiment performed by the National Railway Corporation of Japan.
There were 16 runs and 9 experimental factors. The response variable is the observed
tensile strength of the weld, one of several quality characteristics measured. All other
variables are at plus and minus levels.
}
\usage{data(welding)}
\format{
A data frame containing the variables.
All the explanatory variables are numeric with two levels, \code{-1} and \code{1}.
  \item{Rods}{Kind of welding rods}
  \item{Drying}{Period of drying}
  \item{Material}{Welded material}
  \item{Thickness}{Thickness}
  \item{Angle}{Angle}
  \item{Opening}{Opening}
  \item{Current}{Current}
  \item{Method}{Welding method}
  \item{Preheating}{Preheating}
  \item{Strength}{Tensile strength of the weld in kg/mm. The response variable.}
}
\source{
\url{http://www.statsci.org/data/general/welding.html}
}
\references{
Smyth, G. K., Huele, F., and Verbyla, A. P. (2001). Exact and approximate REML for heteroscedastic regression. \emph{Statistical Modelling} \bold{1}, 161-175.

Smyth, G. K. (2002). An efficient algorithm for REML in heteroscedastic regression. \emph{Journal of Computational and Graphical Statistics} \bold{11}, 1-12.
}
\keyword{datasets}
