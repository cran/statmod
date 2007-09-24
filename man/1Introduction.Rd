\name{1.StatMod}
\alias{1.StatMod}
\title{Introduction to the StatMod Package}
\description{
This library packages together those functions, other than those for microarray data analysis, which I wish to make public.
A change-log for this package is available from
\url{http://www.statsci.org/r/changelog.txt}.

Contributions to this library have also been made by Paul Bagshaw, Centre National d'Etudes des Telecommunications (DIH/DIPS), France (\code{qinvgauss}), and Trevor Park, Department of Statistics, University of Florida (\code{rinvgauss}).
}

\section{Generalized Linear Models}{
\code{tweedie},
\code{canonic.digamma},
\code{unitdeviance.digamma},
\code{varfun.digamma},
\code{cumulant.digamma},
\code{d2cumulant.digamma},
\code{meanval.digamma}
and \code{logmdigamma}
are functions to fit non-standard generalized linear models related to the gamma distribution. 

\code{qres} implements randomized quantile residuals for generalized linear models.
}

\section{Growth Curves}{
\code{compareGrowthCurves},
\code{compareTwoGrowthCurves} and
\code{meanT}
are functions to test for differences between growth curves with repeated measurements on subjects.
}

\section{Limiting Dilution Analysis}{
\code{\link{limdil}}
implements limiting dilution analysis using complemenary log-log binomial generalized linear model regression, with some improvements on previous programs.
}

\section{Probability Distributions}{
\code{qinvgauss},
\code{dinvgauss},
\code{pinvgauss} and
\code{rinvgauss}
perform probability calculations for the inverse Gaussian distribution. 

\code{gauss.quad} and
\code{gauss.quad.prob} compute Gaussian Quadrature with probability distributions. 
}

\section{Tests}{
\code{hommel.test} performs Hommel's multiple comparison tests.

\code{power.fisher.test} computes the power of Fisher's Exact Test for comparing proportions.

\code{sage.test} is a fast approximation to Fisher's exact test for each tag for comparing two Serial Analysis of Gene Expression (SAGE) libraries.
}

\section{Variance Models}{
\code{mixedModel2},
\code{mixedModel2Fit} and
\code{glmgam.fit} fit mixed linear models.

\code{remlscore} and \code{remlscoregamma} fit heteroscedastic and varying dispersion models by REML.
\code{welding} is an example data set.
}

\section{Matrix Computations}{
\code{matvec} and \code{vecmat} facilitate multiplying matrices by vectors.
}

\author{Gordon Smyth}
\keyword{documentation}
