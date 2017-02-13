library(xtable)
library(mlbench)
ozone = data(Ozone,package='mlbench')

ozone = na.omit(Ozone)[,4:13]
y = ozone[,1]
x = as.matrix(ozone[,2:10])
x = cbind(1,x)

#Problem A, Quantifying Uncertainty
betahat = solve(t(x) %*% x) %*% t(x) %*% y
lm1 = lm(y~x-1)
summary(lm1)
 
RSS <- (1/(length(y)-dim(x)[2]))*sum((y - x%*%betahat)^2)
d.matrix <- solve(t(x) %*% x)
sqrt(diag(RSS*d.matrix))

#Problem A, Bootstrap
n.iter <- 10000
beta.mat <- NULL

#xy pairs
for (i in 1:n.iter){
	samp <- sample(1:length(y),length(y),replace=TRUE)

	y.boot    <- y[samp]
	x.boot    <- x[samp,]
	beta.boot <- solve(t(x.boot) %*% x.boot) %*% t(x.boot) %*% y.boot

	beta.mat  <- cbind(beta.mat,beta.boot)
}

cov(t(beta.mat))

#residuals
n.iter <- 10000
beta.mat <- NULL

res     <- y - x%*%betahat
y.hat   <- y-res
hat.mat <- solve(t(x) %*% x) %*% t(x) 

for (i in 1:n.iter){
	samp <- sample(res,length(y),replace=TRUE)

	y.boot    <- y.hat - samp
	beta.boot <- hat.mat%*% y.boot

	beta.mat  <- cbind(beta.mat,beta.boot)
}

cov(t(beta.mat))

#Covariance matrix from asymptotic theory (Calculated previously)
RSS*d.matrix

#Problem B, Bootstrap

A <- matrix(c(4,2,2,3),2,2)
b <- matrix(c(-1,1),2,1)
w <- matrix(0,2,500)

v<-eigen(A)

# Draw 500 bivariate normals
for (j in 1:500){
	x <- matrix(rnorm(2),2,1) 
	draws <- v$vectors%*%diag(sqrt(v$values))%*%solve(v$vectors)%*%x+b
	w[1,j] <- draws[1]
	w[2,j] <- draws[2]
}

# Estimate mean and variance
rowSums(w)/500
#[1] -1.0900930  0.9165089
cov(t(w))
#       [,1]     [,2]
#[1,] 4.044331 1.922020
#[2,] 1.922020 2.830794

# Bootstrap sample
dim.mat <- matrix(0,n.iter,2)
for (i in 1:n.iter){
	samp <- sample(1:500,500,replace=TRUE)

	w1    <- w[1,samp]
	w2    <- w[2,samp]

	dim.mat[i,1]  <- sum(w1)/500
	dim.mat[i,2]  <- sum(w2)/500
}

cov(dim.mat)
#            [,1]        [,2]
#[1,] 0.007823330 0.003771656
#[2,] 0.003771656 0.005536346