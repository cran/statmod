library(statmod)

set.seed(0); u <- runif(100)

### glmgam.fit

glmgam.fit(1,1)
glmgam.fit(c(1,1),c(0,4))
glmgam.fit(X=cbind(1,c(1,0.5,0.5,0,0)),y=rchisq(5,df=1))

### randomizedBlock

y <- rnorm(6)
x <- rnorm(6)
z <- c(1,1,2,2,3,3)
randomizedBlock(y~x,random=z)

### randomizedBlockFit

y <- c(-1,1,-2,2,0.5)
X <- matrix(1,5,1)
Z <- model.matrix(~factor(c(1,1,2,2,3))-1)
randomizedBlockFit(y,X,Z)

### qresiduals

y <- rnorm(6)
fit <- glm(y~1)
residuals(fit)
qresiduals(fit)
qresiduals(fit,dispersion=1)
