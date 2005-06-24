#  LIMDIL.R

limdil <- function(response,dose,tested=rep(1,length(response)),observed=FALSE,confidence=0.95,test.unit.slope=FALSE)
#	Limiting dilution analysis
#	Gordon Smyth
#	21 June 2005
{
	alpha <- 1-confidence
	out <- list()
	f <- binomial(link="cloglog")
	f$aic <- quasi()$aic
	y <- response/tested
	if(any(y < 0)) stop("Negative values for response or tested")
	if(any(y > 1)) stop("response > tested")
	if(all(y < 1e-15)) {
		N <- sum(dose*tested)
		if(observed)
			U <- 1-alpha^(1/N)
		else
			U <- -log(alpha)/N
		out$CI <- c(Lower=Inf,Estimate=Inf,Upper=1/U)
		if(test.unit.slope) {
			out$test.unit.slope <- c(Chisq=NA,P.value=NA)
		}
		return(out)
	}
	if(all(1-y < 1e-15)) {
		dosem <- min(dose)
		tested <- sum(tested[dose==dosem])
		dose <- dosem
		beta <- log(-log(1-alpha^(1/tested)))-log(dose)
		if(observed)
			U <- 1-exp(-exp(beta))
		else
			U <- exp(beta)
		out$CI <- c(Lower=1/U,Estimate=0,Upper=0)
		if(test.unit.slope) {
			out$test.unit.slope <- c(Chisq=NA,P.value=NA)
		}
		return(out)
	}
	fit0 <- glm(y~offset(log(dose)),family=f,weights=tested)
	s <- summary(fit0)
	Estimate <- s$coef[,"Estimate"]
	SE <- s$coef[,"Std. Error"]
	z <- qnorm(alpha/2,lower.tail=FALSE)
	CI <- c(Lower=Estimate-z*SE,Estimate=Estimate,Upper=Estimate+z*SE)
	if(observed)
		Frequency <- 1-exp(-exp(CI))
	else
		Frequency <- exp(CI)
	out$CI <- 1/Frequency
	if(test.unit.slope && length(response)>1) {
		fit1 <- glm(y~log(dose),family=f,weights=tested)
		dev <- fit0$deviance-fit1$deviance
		p.value <- pchisq(dev,df=1,lower.tail=FALSE)
		out$test.unit.slope <- c(Chisq=dev,P.value=p.value)
	}
	out
}

