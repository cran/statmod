#  SAGE.R

fisher.test2 <- function(x, workspace=2e+05)
#	Sped up version of ctest: fisher.test
#	Gordon Smyth
#	15 Nov 2003
{
	if(!is.matrix(x)) stop("x must be a matrix")
	d <- dim(x)
	if(any(d < 2)) stop("x must have at least 2 rows and columns")
	if(any(x < 0) || any(is.na(x))) stop("all entries of x must be nonnegative and finite")
	nr <- d[1]
	nc <- d[2]
	.C("fexact", as.integer(nr), as.integer(nc), 
		as.double(x), as.integer(nr), as.double(-1), as.double(100), 
		as.double(0), as.double(0), p = as.double(0), as.integer(workspace), 
		PACKAGE = "ctest")$p
}

sage.test <- function(x, y, n1=sum(x), n2=sum(y))
#	Binomial probabilities for comparing SAGE libraries
#	Gordon Smyth
#	15 Nov 2003.  Last modified 18 Nov 2003.
{
	if(any(is.na(x)) || any(is.na(y))) stop("missing values not allowed")
	x <- as.integer(x)
	y <- as.integer(y)
	if(any(x<0) || any(y<0)) stop("x and y must be non-negative")
	if(length(x) != length(y)) stop("x and y must have same length")
	size <- x+y
	if(n1==n2) {
		x <- pmin(x,y)
		return(pbinom(x,size=size,prob=0.5)+pbinom(size-x+0.5,size=size,prob=0.5,lower.tail=FALSE))
	}
	prob <- n1/(n1+n2)
	p.value <- rep(1,length(x))
	if(any(big <- size>10000)) {
		ibig <- (1:length(x))[big]
		for (i in ibig) p.value[i] <- chisq.test(matrix(c(x[i],y[i],n1-x[i],n2-y[i]),2,2))$p.value
	}
	size0 <- size[size>0 & !big]
	if(length(size0)) for (isize in unique(size0)) {
		i <- (size==isize)
		p <- dbinom(0:isize,p=prob,size=isize)
		o <- order(p)
		cumsump <- cumsum(p[o])[order(o)]
		p.value[i] <- cumsump[x[i]+1]
	}
	p.value
}
