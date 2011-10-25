library(statmod)

set.seed(0); u <- runif(100)

### fitNBP

y <- matrix(rnbinom(2*4,mu=4,size=1.5),2,4)
lib.size <- rep(50000,4)
group <- c(1,1,2,2)
fitNBP(y,group=group,lib.size=lib.size)

### glmgam.fit

glmgam.fit(1,1)
glmgam.fit(c(1,1),c(0,4))
glmgam.fit(X=cbind(1,c(1,0.5,0.5,0,0)),y=rchisq(5,df=1))

### glmnb.fit
y <- rnbinom(5,mu=10,size=10)
glmnb.fit(X=cbind(1,c(1,0.5,0.5,0,0)),y=y,dispersion=0.1)
glmnb.fit(X=cbind(1,c(1,0.5,0.5,0,0)),y=y,dispersion=runif(6))
glmnb.fit(X=cbind(1,c(1,1,0,0,0)),y=c(0,0,6,2,9),dispersion=0.1)

### mixedModel2

y <- rnorm(6)
x <- rnorm(6)
z <- c(1,1,2,2,3,3)
mixedModel2(y~x,random=z)

### mixedModel2Fit

y <- c(-1,1,-2,2,0.5)
X <- matrix(1,5,1)
Z <- model.matrix(~factor(c(1,1,2,2,3))-1)
mixedModel2Fit(y,X,Z)

### qresiduals

y <- rnorm(6)
fit <- glm(y~1)
residuals(fit)
qresiduals(fit)
qresiduals(fit,dispersion=1)

if(require("MASS")) {
fit <- glm(Days~Age,family=negative.binomial(2),data=quine)
print(summary(qresiduals(fit)))
fit <- glm.nb(Days~Age,link=log,data = quine)
print(summary(qresiduals(fit)))
}

