#  NUMERICAL INTEGRATION

gauss.quad <- function(n,kind="legendre",alpha=0,beta=0) {
#	Calculate nodes and weights for Guassian quadrature.
#	Adapted from Netlib routine gaussq.f
#	Gordon Smyth, Walter and Eliza Hall Institute, smyth@wehi.edu.au
#	4 Sept 2002

	kind <- match.arg(kind,c("legendre","chebyshev1","chebyshev2","hermite","jacobi","laguerre"))
	i <- 1:n
	i1 <- 1:(n-1)
	switch(kind, legendre={
		muzero <- 2
		a <- rep(0,n)
		b <- i1/sqrt(4*i1^2-1)
	}, chebyshev1={
		muzero <- pi
		a <- rep(0,n)
		b <- rep(0.5,n-1)
		b[1] <- sqrt(0.5)
	}, chebyshev2={
		muzero <- pi/2
		a <- rep(0,n)
		b <- rep(0.5,n-1)
	}, hermite={
		muzero <- sqrt(pi)
		a <- rep(0,n)
		b <- sqrt(i1/2)
	}, jacobi={
		ab <- alpha+beta
		muzero <- 2^(ab+1)*gamma(alpha+1)*gamma(beta+1)/gamma(ab+2)
		a <- i
		a[1] <- (beta-alpha)/(ab+2)
		i2 <- 2:n
		abi <- ab+2*i2
		a[i2] <- (beta^2-alpha^2)/(abi-2)/abi
		b <- i1
		b[1] <- sqrt(4*(alpha+1)*(beta+1)/(ab+2)^2/(ab+3))
		i2 <- 2:(n-1)
		abi <- ab+2*i2
		b[i2] <- sqrt(4*i2*(i2+alpha)*(i2+beta)*(i2+ab)/(abi^2-1)/abi^2)
	}, laguerre={
		a <- 2*i-1+alpha
		b <- sqrt(i1*(i1+alpha))
		muzero <- gamma(alpha+1)
	})
	A <- rep(0,n*n)
	A[(n+1)*(i-1)+1] <- a
	A[(n+1)*(i1-1)+2] <- b
	A[(n+1)*i1] <- b
	dim(A) <- c(n,n)
	vd <- eigen(A,symmetric=TRUE)
	w <- rev(as.vector( vd$vectors[1,] ))
	w <- muzero * w^2
	x <- rev( vd$values )
	list(nodes=x,weights=w)
}

gauss.quad.prob <- function(n,dist="uniform",l=0,u=1,mu=0,sigma=1,alpha=1,beta=1) {
#	Calculate nodes and weights for Guassian quadrature using probability densities.
#	Adapted from Netlib routine gaussq.f
#	Gordon Smyth, Walter and Eliza Hall Institute, smyth@wehi.edu.au
#	4 Sept 2002

	dist <- match.arg(dist,c("uniform","beta1","beta2","normal","beta","gamma"))
	if(dist=="beta" && alpha==0.5 && beta==0.5) dist <- "beta1"
	if(dist=="beta" && alpha==1.5 && beta==1.5) dist <- "beta2"
	i <- 1:n
	i1 <- 1:(n-1)
	switch(dist, uniform={
		a <- rep(0,n)
		b <- i1/sqrt(4*i1^2-1)
	}, beta1={
		a <- rep(0,n)
		b <- rep(0.5,n-1)
		b[1] <- sqrt(0.5)
	}, beta2={
		a <- rep(0,n)
		b <- rep(0.5,n-1)
	}, normal={
		a <- rep(0,n)
		b <- sqrt(i1/2)
	}, beta={
		ab <- alpha+beta
		a <- i
		a[1] <- (alpha-beta)/ab
		i2 <- 2:n
		abi <- ab-2+2*i2
		a[i2] <- ((alpha-1)^2-(beta-1)^2)/(abi-2)/abi
		b <- i1
		b[1] <- sqrt(4*alpha*beta/ab^2/(ab+1))
		i2 <- 2:(n-1)
		abi <- ab-2+2*i2
		b[i2] <- sqrt(4*i2*(i2+alpha-1)*(i2+beta-1)*(i2+ab-2)/(abi^2-1)/abi^2)
	}, gamma={
		a <- 2*i+alpha-2
		b <- sqrt(i1*(i1+alpha-1))
	})
	A <- rep(0,n*n)
	A[(n+1)*(i-1)+1] <- a
	A[(n+1)*(i1-1)+2] <- b
	A[(n+1)*i1] <- b
	dim(A) <- c(n,n)
	vd <- eigen(A,symmetric=TRUE)
	w <- rev(as.vector( vd$vectors[1,] ))^2
	x <- rev( vd$values )
	switch(dist,
		uniform = x <- l+(u-l)*(x+1)/2,
		beta1 = x <- (x+1)/2,
		beta2 = x <- (x+1)/2,
		normal = x <- mu + sqrt(2)*sigma*x,
		beta = x <- (x+1)/2,
		gamma = x <- beta*x)
	list(nodes=x,weights=w)
}
