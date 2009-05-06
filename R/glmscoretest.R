##  glmscore.R

glm.scoretest <- function(fit, x2)
#	Score test for new covariate in glm
#	Gordon Smyth
#	27 March 2009
{
	ws <- sqrt(fit$weights)
	x2.1w <- qr.resid(fit$qr,ws*x2)
	zw <- ws*fit$residuals
	colSums(as.matrix(x2.1w*zw))/sqrt(colSums(as.matrix(x2.1w * x2.1w)))
}

