library(MASS)

# Squared Exponential

sqexp.cov <- function(tau1,tau2,b,x){
	temp <- rep(x,length(x))
	dist <- abs(matrix(temp,length(x),length(x))-matrix(temp,length(x),length(x),byrow=TRUE))	
	cov  <- (tau1^2)*exp(-(1/2)*(dist/b)^2)+diag(rep(tau2^2,length(x)))
	cov
} 

# Matern

matern52.cov <- function(tau1,tau2,b,x){
	temp <- rep(x,length(x))
	dist <- abs(matrix(temp,length(x),length(x))-matrix(temp,length(x),length(x),byrow=TRUE))	
	cov  <- (tau1^2)*(1 + sqrt(5)*dist/b + 5*(dist^2)/(3*b^2))*exp(-(sqrt(5)*dist/b)) + diag(rep(tau2^2,length(x)))
	cov	
}


X <- runif(100)

par(mfrow=c(3,3))

# Constant Tau2

# Row

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.0001,0.000001,0.001,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.0001,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.0001,0.000001,0.1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.0001,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.0001,0.000001,1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.0001,0.000001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.001,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.01,0.000001,0.001,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.01,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.01,0.000001,0.1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.01,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.01,0.000001,1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.01,0.000001,1,X)),col=i,pch=20)}


# Constant Tau1

# Row

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.001,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.000001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.001,0.001,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.001,0.1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.001,1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.01,0.001,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.01,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.01,0.1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.01,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.01,1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),sqexp.cov(0.001,0.01,1,X)),col=i,pch=20)}

############## Matern ################

# Constant Tau2

# Row

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.0001,0.000001,0.001,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.0001,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.0001,0.000001,0.1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.0001,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.0001,0.000001,1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.0001,0.000001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.001,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.01,0.000001,0.001,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.01,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.01,0.000001,0.1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.01,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.01,0.000001,1,X)),ylab='Y',ylim=c(-0.03,0.03),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.01,0.000001,1,X)),col=i,pch=20)}


# Constant Tau1

# Row

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.001,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.000001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.001,0.001,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.001,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.001,0.1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.001,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.001,1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.001,1,X)),col=i,pch=20)}

# Row

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.01,0.001,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.01,0.001,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.01,0.1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.01,0.1,X)),col=i,pch=20)}

plot(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.01,1,X)),ylab='Y',ylim=c(-0.02,0.02),pch=20)
for (i in 2:4){
points(X,mvrnorm(1,rep(0,length(X)),matern52.cov(0.001,0.01,1,X)),col=i,pch=20)}


