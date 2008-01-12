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

### limdil

#####################################
# Test 1, when there is only 1 group
#####################################
Dose <- c(50,100,200,400,800)
Responses <- c(2,6,9,15,21)
Tested <- c(24,24,24,24,24)
Group <- c(1,1,1,1,1)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)

###################################
# Test 2, when there is 1 group 
# 1 with no responses
###################################
Dose <- c(30000,20000,4000,500)
Responses <- c(0,0,0,0)
Tested <- c(6,6,6,6)
Group <- c(1,1,1,1)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)

###################################
# Test 3, when there is 1 group 
# 1 with all responses
###################################
Dose <- c(30000,20000,4000,500)
Responses <- c(6,6,6,6)
Tested <- c(6,6,6,6)
Group <- c(1,1,1,1)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)

##################################
# Test 4, when there are 2 groups
# average response
##################################
Dose <- c(30000,20000,4000,500,30000,20000,4000,500)
Responses <- c(2,3,2,3,2,2,4,3)
Tested <- c(6,6,6,6,6,6,6,6)
Group <- c(1,1,1,1,2,2,2,2)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)

##################################
# Test 5, when there are 2 groups
# 1 with nearly all responses
# 1 with nearly no responses
##################################
Dose <- c(30000,20000,4000,500,30000,20000,4000,500)
Responses <- c(1,0,0,0,6,6,6,5)
Tested <- c(6,6,6,6,6,6,6,6)
Group <- c(1,1,1,1,2,2,2,2)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)

###################################
# Test 6, when there is 2 groups, 
# 1 with all responses
# 1 with none responses
###################################
Dose <- c(30000,20000,4000,500,30000,20000,4000,500)
Responses <- c(0,0,0,0,6,6,6,6)
Tested <- c(6,6,6,6,6,6,6,6)
Group <- c(1,1,1,1,2,2,2,2)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)

##################################
# Test 7, when there are 4 groups
# 4 average response
##################################
Dose <- c(30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500)
Responses <- c(2,3,2,1,6,5,6,1,2,3,4,2,6,6,6,1)
Tested <- c(6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6)
Group <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)

##################################
# Test 8, when there are 4 groups 
# 1 group with no response
# 1 group with all response
# 2 group average response
##################################
Dose <- c(30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500,30000,20000,4000,500)
Responses <- c(0,0,0,0,6,5,6,1,2,3,4,2,6,6,6,6)
Tested <- c(6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6)
Group <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)

limdil(Responses,Dose,Tested,Group,test.unit.slope=TRUE)
