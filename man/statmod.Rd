\name{statmod-package}
\alias{statmod}
\alias{statmod-package}
\docType{package}
\title{Introduction to the StatMod Package}
\description{
This package includes a variety of functions for numerical analysis and statistical modelling.
The functions are briefly summarized by type of application below.
}

\section{Generalized Linear Models}{
The function \code{\link{tweedie}} defines a large class of generalized linear model families with power variance functions.
It used in conjunction with the glm function, and widens the class of families that can be fitted.

\code{\link{qresiduals}} implements randomized quantile residuals for generalized linear models.

The functions
\code{canonic.digamma},
\code{unitdeviance.digamma},
\code{varfun.digamma},
\code{cumulant.digamma},
\code{d2cumulant.digamma},
\code{meanval.digamma}
and \code{logmdigamma}
are used to fit double-generalized models, in which a link-linear model is fitted to the dispersion as well as to the mean.
Spefically they are used to fit the dispersion submodel associated with a gamma glm.
}

\section{Growth Curves}{
\code{\link{compareGrowthCurves}},
\code{compareTwoGrowthCurves} and
\code{meanT}
are functions to test for differences between growth curves with repeated measurements on subjects.
}

\section{Limiting Dilution Analysis}{
The \code{\link{limdil}} function is used in the analysis of stem cell frequencies.
It implements limiting dilution analysis using complemenary log-log binomial generalized linear model regression, with some improvements on previous programs.
}

\section{Probability Distributions}{
The functions
\code{\link{qinvgauss}},
\code{\link{dinvgauss}},
\code{\link{pinvgauss}} and
\code{\link{rinvgauss}}
provide probability calculations for the inverse Gaussian distribution. 

\code{\link{gauss.quad}} and
\code{gauss.quad.prob} compute Gaussian Quadrature with probability distributions. 
}

\section{Tests}{
\code{\link{hommel.test}} performs Hommel's multiple comparison tests.

\code{\link{power.fisher.test}} computes the power of Fisher's Exact Test for comparing proportions.

\code{\link{sage.test}} is a fast approximation to Fisher's exact test for each tag for comparing two Serial Analysis of Gene Expression (SAGE) libraries.

\code{\link{permp}} computes p-values for permutation tests when the permutations are randomly drawn.
}

\section{Variance Models}{
\code{\link{mixedModel2}},
\code{\link{mixedModel2Fit}} and
\code{\link{glmgam.fit}} fit mixed linear models.

\code{\link{remlscore}} and \code{\link{remlscoregamma}} fit heteroscedastic and varying dispersion models by REML.
\code{\link{welding}} is an example data set.
}

\section{Matrix Computations}{
\code{\link{matvec}} and \code{\link{vecmat}} facilitate multiplying matrices by vectors.
}

\author{Gordon Smyth}
\keyword{documentation}
