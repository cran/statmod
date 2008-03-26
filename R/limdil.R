#  LIMDIL.R

limdil <- function (response, dose, tested = rep(1, length(response)), group=rep(1,length(response)),observed = FALSE, confidence = 0.95, test.unit.slope = FALSE) 
#	Limiting dilution analysis
#	Gordon Smyth, Yifang Hu
#	21 June 2005. Last revised 05 February 2008
{
#	Check input arguments
	group <- as.factor(group)

	alpha <- 1 - confidence
	out <- list()
	f <- binomial(link = "cloglog")
	f$aic <- quasi()$aic
	y <- response/tested
	if (any(y < 0)) stop("Negative values for response or tested")
	if (any(y > 1)) stop("The response cannot be greater than the number tested")

	num.group<-length(levels(group))

	if(num.group==1)
	{
		fit0 <- glm(y ~ offset(log(dose)), family = f, weights = tested)
		s <- summary(fit0)
		Estimate <- s$coef[, "Estimate"]
		SE <- s$coef[, "Std. Error"]
		z <- qnorm(alpha/2, lower.tail = FALSE)
		CI <- c(Lower = Estimate - z * SE, Estimate = Estimate, Upper = Estimate + z * SE)

		if (observed) Frequency <- 1 - exp(-exp(CI))
		else Frequency <- exp(CI)
		out$CI <- 1/Frequency
		if (test.unit.slope && length(response) > 1) 
		{
		fit1 <- glm(y ~ log(dose), family = f, weights = tested)
		dev <- fit0$deviance - fit1$deviance
		p.value <- pchisq(dev, df = 1, lower.tail = FALSE)
		out$test.unit.slope <- c(Chisq = dev, P.value = p.value, df=1)
		}
	}

	if(num.group>1)
	{
		fit0 <- glm(y ~ offset(log(dose))+group-1, family = f, weights = tested)
		s <- summary(fit0)
		Estimate <- s$coef[, "Estimate"]
		SE <- s$coef[, "Std. Error"]
		z <- qnorm(alpha/2, lower.tail = FALSE)
		CI <- cbind(Lower = Estimate - z * SE, Estimate = Estimate, Upper = Estimate + 
		z * SE)
		if (observed) Frequency <- 1 - exp(-exp(CI))
		else Frequency <- exp(CI)
		out$CI <- 1/Frequency
		rownames(out$CI)<-paste("Group",levels(group))

		fit2 <- glm(y~offset(log(dose))+group,family=f,weights=tested)
		dev.g<-fit2$null.deviance-fit2$deviance
		group.p<-pchisq(dev.g,df=num.group-1,lower=FALSE)
		out$test.difference<-c(Chisq = dev.g, P.value = group.p, df=num.group-1)

		if (test.unit.slope && length(response) > 1) 
		{
			fit3 <- glm(y ~ log(dose)+group, family = f, weights = tested)
			dev <- fit2$deviance - fit3$deviance
			slope.p <- pchisq(dev, df = 1, lower.tail = FALSE)
			out$test.unit.slope <- c(Chisq = dev, P.value = slope.p, df = 1)
		}
	}
	
	groupLevel<-levels(group)

	for(i in 1:num.group)
	{
	 	tmp <- (group == groupLevel[i])
		index<-which(tmp==TRUE)
		y.element<-response[index]/tested[index]	

		if(all(y.element<1e-15))
		{
		 	N <- sum(dose * tested)
			if (observed) U <- 1 - alpha^(1/N)
			else U <- -log(alpha)/N
			
			if(num.group==1)
			{
				out$CI<-c(Lower = Inf, Estimate = Inf, Upper = 1/U)
				if (test.unit.slope)
				{
					out$test.unit.slope<- c(Chisq = NA, P.value = NA, df=1)
				}
			}

			else out$CI[i,]<-c(Lower = Inf, Estimate = Inf, Upper = 1/U)
		}
		
		if(all(1- y.element<1e-15))
		{
			dosem <- min(dose)
			tested.sum <- sum(tested[dose == dosem])
			dose <- dosem
			beta <- log(-log(1 - alpha^(1/tested.sum))) - log(dose)
			if (observed) U <- 1 - exp(-exp(beta))
			else U <- exp(beta)
			
			if(num.group==1)
			{
				out$CI<- c(Lower = 1/U, Estimate = 0, Upper = 0)
				if (test.unit.slope)
				{
						out$test.unit.slope<- c(Chisq = NA, P.value = NA, df=1)
				}

			}
			else out$CI[i,]<- c(Lower = 1/U, Estimate = 0, Upper = 0)
		}
	}
	return(out)
}