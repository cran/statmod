glmgam.fit <- function(X,y,start=NULL,trace=FALSE,tol=1e-6,maxit=50) {
#  Fit gamma generalized linear model with identity link
#  by Levenberg damped Fisher scoring
#  Gordon Smyth, smyth@wehi.edu.au
#  12 Mar 2003.  Last revised 13 Mar 2003.

#  check input
X <- as.matrix(X)
y <- as.vector(y)
n <- length(y)
if(n != nrow(X)) stop("length(y) not equal to nrow(X)")
if(n == 0) return(list(coefficients=numeric(0),fitted.values=numeric(0),deviance=numeric(0)))
if(any(is.na(y)) || any(is.na(X))) stop("Missing values not allowed")
if(any(y <= 0)) stop("y must be positive")
p <- ncol(X)
if(p > n) stop("More columns than rows in X")

#  starting values
if(is.null(start)) {
	beta <- lm.wfit(X,y,1/y^2)$coef
	mu <- X %*% beta
	if(any(mu <= 0)) {
		beta <- lm.fit(X,rep(mean(y),n))$coef
		mu <- X %*% beta
		if(any(mu <= 0)) stop("Can't find admissable starting values")
	}
} else {
	beta <- start
	mu <- X %*% beta
	if(any(mu <= 0)) stop("Starting values not admissable")
}
dev <- 2*sum( (y-mu)/mu - log(y/mu) )

# reml scoring
iter <- 0
if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")
repeat {
	iter <- iter+1

	# information matrix
	v <- mu^2
	XVX <- crossprod(X,vecmat(1/v,X))
	maxinfo <- max(diag(XVX))
	if(iter==1) {
		lambda <- abs(mean(diag(XVX)))/p
		I <- diag(p)
	}

	# score vector
	dl <- crossprod(X,(y-mu)/v)

	# Levenberg damping
	betaold <- beta
	devold <- dev
	lev <- 0
	repeat {
		lev <- lev+1

		# trial step
		R <- chol(XVX + lambda*I)
		dbeta <- backsolve(R,backsolve(R,dl,transpose=TRUE))
		beta <- betaold + dbeta
		mu <- X %*% beta
		if(all(mu>0))
			dev <- 2*sum( (y-mu)/mu - log(y/mu) )
		else
			dev <- Inf
		if(is.finite(dev)) if(dev < devold || dev/max(mu) < 1e-15) break

		# exit if too much damping
		if(lambda/maxinfo > 1e15) {
			beta <- betaold
			warning("Too much damping - convergence tolerance not achievable")
			break
		}

		# step not successful so increase damping
		lambda <- 2*lambda
		if(trace) cat("Damping increased to",lambda,"\n")
	}

	# iteration output
	if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")

	# keep exiting if too much damping
	if(lambda/maxinfo > 1e15) break

	# decrease damping if successful at first try
	if(lev==1) lambda <- lambda/10

	# test for convergence
	if( crossprod(dl,dbeta) < tol || dev/max(mu) < 1e-15) break

	# test for iteration limit
	if(iter > maxit) {
		warning("Max iterations exceeded")
		break
	}
}

list(coefficients=as.vector(beta),fitted.values=as.vector(mu),deviance=dev)
}
