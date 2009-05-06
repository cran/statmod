#  LIMDIL.R

elda <- limdil <- function (response, dose, tested = rep(1, length(response)), group=rep(1,length(response)),observed = FALSE, confidence = 0.95, test.unit.slope = FALSE) 
#	Limiting dilution analysis
#	Gordon Smyth, Yifang Hu
#	21 June 2005. Last revised 27 March 2009.
{
	group <- as.factor(group)

	alpha <- 1 - confidence
	out <- list()
	f <- binomial(link = "cloglog")
	f$aic <- quasi()$aic
	y <- response/tested
	if (any(y < 0)) stop("Negative values for response or tested")
	if (any(y > 1)) stop("The response cannot be greater than the number tested")

	num.group<-length(levels(group))
	groupLevel<-levels(group)

	out$response <- response
	out$tested <- tested
	out$dose <- dose
	out$group <- group
	out$num.group <- num.group
	class(out) <- "limdil"
	
	if(num.group==1)
	{
		tmp <- (group == groupLevel[1])
		index<-which(tmp==TRUE)
		y.element<-response[index]/tested[index]

		if(all(y.element<1e-15))
		{
		 	N <- sum(dose[index] * tested[index])
			if (observed) U <- 1 - alpha^(1/N)
			else U <- -log(alpha)/N
			out$CI<-c(Lower = Inf, Estimate = Inf, Upper = 1/U)
			if (test.unit.slope) out$test.unit.slope<- c(Chisq = NA, P.value = NA, df=1)
	   
			return(out)
		}
		
		if(all(1- y.element<1e-15))
		{
			U <- .limdil.allpos(tested=tested[index],dose=dose[index],confidence=confidence,observed=observed)
			out$CI<- c(Lower = 1/U, Estimate = 1, Upper = 1)
			if (test.unit.slope) out$test.unit.slope<- c(Chisq = NA, P.value = NA, df=1)

			return(out)
		}
	
	
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
		y.element.all<-response/tested
		
		if(all(y.element.all<1e-15) || all(1-y.element.all<1e-15))
		{
			out$CI<-matrix(nrow=num.group,ncol=3)
			colnames(out$CI)<-c("Lower","Estimate","Upper")
			rownames(out$CI)<-paste("Group",groupLevel)
			
			for(i in 1:num.group)
			{
	 			tmp <- (group == groupLevel[i])
				index<-which(tmp==TRUE)
				y.element<-response[index]/tested[index]
				
				if(all(y.element<1e-15))
				{
					N <- sum(dose[index] * tested[index])
					if (observed) U <- 1 - alpha^(1/N)
					else U <- -log(alpha)/N
					out$CI[i,]<-c(Lower = Inf, Estimate = Inf, Upper = 1/U)
				}
					
				if(all(1- y.element<1e-15))
				{
					U <- .limdil.allpos(tested=tested[index],dose=dose[index],confidence=confidence,observed=observed)
					out$CI[i,]<- c(Lower = 1/U, Estimate = 1, Upper = 1)
					
				}
			}	
			out$test.difference<-c(Chisq = 0, P.value = 1, df=num.group-1)
			if (test.unit.slope) out$test.unit.slope<- c(Chisq = NA, P.value = NA, df=1)
			
			return(out)
		}
	
		fit0 <- glm(y ~ offset(log(dose))+group-1, family = f, weights = tested)
		s <- summary(fit0)
		Estimate <- s$coef[, "Estimate"]
		SE <- s$coef[, "Std. Error"]
		z <- qnorm(alpha/2, lower.tail = FALSE)
		CI <- cbind(Lower = Estimate - z * SE, Estimate = Estimate, Upper = Estimate + z * SE)
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
	

	for(i in 1:num.group)
	{
	 	tmp <- (group == groupLevel[i])
		index<-which(tmp==TRUE)
		y.element<-response[index]/tested[index]

		if(all(y.element<1e-15))
		{
		 	N <- sum(dose[index] * tested[index])
			if (observed) U <- 1 - alpha^(1/N)
			else U <- -log(alpha)/N
			out$CI[i,]<-c(Lower = Inf, Estimate = Inf, Upper = 1/U)
		}

		if(all(1- y.element<1e-15))
		{
			U <- .limdil.allpos(tested=tested[index],dose=dose[index],confidence=confidence,observed=observed)
			out$CI[i,]<- c(Lower = 1/U, Estimate = 1, Upper = 1)
		}
	}

	return(out)

}

print.limdil <- function(x, ...)
{
#	Print limiting dilution analysis
#	Yifang Hu and Gordon Smyth
#	20 February 2009. Last revised 27 March 2009.

	cat("Confidence intervals for frequency:\n\n")
	print(x$CI)
	cat("\n\n")
	difference<-NULL
	single<-NULL

	if(is.null(x$test.difference)!=TRUE) difference<-x$test.difference
	if(is.null(x$test.unit.slope)!=TRUE) single<-x$test.unit.slope

	if((is.null(x$test.difference)!=TRUE) || (is.null(x$test.unit.slope)!=TRUE))
	{
		a <- data.frame(rbind(difference,single))
		a <- a[,c(1,3,2)]
		colnames(a) <- c("Chisq","Df","Pr(>Chisq)")

		if(is.null(difference)!=TRUE)
		{	
			rownames(a)[1] <- "Differences between groups"
			if(is.null(single)!=TRUE) rownames(a)[2] <-"Single vs multi-hit"
		}
		else{ rownames(a)<-"Single vs multi-hit"}

		suppressWarnings(printCoefmat(a,tst.ind=1,has.Pvalue=TRUE,P.value=TRUE))
	}
}

plot.limdil <- function(x, col.group=NULL, ...)
#	Plot limiting dilution analysis
#	Yifang Hu and Gordon Smyth
#	20 February 2009. Last revised 23 February 2009.
{
	num.group<-length(levels(x$group))
	if(is.null(col.group)) col.group <- 1:num.group

	col <- x$group
	levels(col) <- col.group
	col <- as.vector(col)
	dose<-x$dose
	maxx<-max(dose)	
	
	if(any(x$response==x$tested)==TRUE)
	{
		i<-which(x$response==x$tested)
		x$response[i]<-x$response[i]-0.5

		nonres<-log(1-x$response/x$tested)
		if(num.group>1) nonres<-pmin(0,jitter(nonres))
		
		miny<-min(nonres)

		plot(dose[-i],nonres[-i],xlim=c(0,maxx),ylim=c(min(miny,-0.5),0),xlab="dose (number of cells)",ylab="log fraction nonresponding", col=col[-i],...)
		points(dose[i],nonres[i],pch=6,col=col[i],...)
	}
	else
	{
		nonres<-log(1-x$response/x$tested)
		miny<-min(nonres)
		if(num.group>1) nonres<-pmin(0,jitter(nonres))		
		plot(dose,nonres,,xlim=c(0,maxx),ylim=c(min(miny,-0.5),0),xlab="dose (number of cells)",ylab="log fraction nonresponding", col=col, ...)
	}

	if(num.group==1)
	{	
		abline(a=0,b=-exp(log(1/(x$CI[2]))),col=col.group, lty=1,...)
		abline(a=0,b=-exp(log(1/(x$CI[1]))),col=col.group, lty=2,...)
		abline(a=0,b=-exp(log(1/(x$CI[3]))),col=col.group, lty=2,...)
	}
	else
	{
		for(i in 1:num.group)
		{
			abline(a=0,b=-exp(log(1/(x$CI[i,2]))),col=col.group[i],lty=1,...)
			abline(a=0,b=-exp(log(1/(x$CI[i,1]))),col=col.group[i],lty=2,...)
			abline(a=0,b=-exp(log(1/(x$CI[i,3]))),col=col.group[i],lty=2,...)

		}
		legend("bottomleft",legend=paste("Group",levels(x$group)),col=col.group,text.col=col.group,cex=0.6, ...)
	}
}

.limdil.allpos <- function(tested, dose, confidence, observed)
#	Yifang Hu. 18 March 2009.
{
	alpha <- 1 - confidence

	dosem <- min(dose)
	tested.group<-tested
	tested.sum <- sum(tested.group[dose == dosem])
	beta <- log(-log(1 - alpha^(1/tested.sum))) - log(dosem)

	if (observed) U <- 1 - exp(-exp(beta))
	else U <- exp(beta)

	lambda <- U
	
	repeat
	{
		if(observed) f <- sum(tested*log(1-(1-lambda)^dose))-log(alpha)
		else f <- sum(tested*log(1-exp(-lambda*dose)))-log(alpha)
		if(observed) deriv <- sum(tested*(-dose)*(1-lambda)^(dose-1)/(1-(1-lambda)^dose)) 
		else deriv <- sum(tested*dose*exp(-dose*lambda)/(1-exp(-dose*lambda)))
		step <- f/deriv
		lambda <- lambda-step
		if(-step < 1e-6)
			break
	}
	lambda
}