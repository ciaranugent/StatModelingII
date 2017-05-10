library(MASS)
gdpgrowth <- read.csv('/Users/ciaranugent/Documents/gdpgrowth.csv')

# Bayesian Linear Model

X <- cbind(1,gdpgrowth$DEF60)
Y <- gdpgrowth$GR6096

K <- diag(.11,2)
m <- matrix(0,2,1)
d <- 1
eta <- 1
omega0 <- 1
mu.star <- t((t(Y)%*%X)%*%solve(t(X)%*%X + K))
lambda.star <- t(X)%*%X + K
eta.star <- eta + sum(Y^2)-t(Y)%*%X%*%solve(t(X)%*%X + K)%*%t(t(Y)%*%X)

n.iter <- 5000

beta <- matrix(0,n.iter,2)
omega <- rep(0,n.iter)

for (i in 1:n.iter){
	beta0 <- mvrnorm(1,mu.star,solve(omega0*lambda.star))
	omega0 <- rgamma(1,(d+dim(X)[1])/2,eta.star/2)
	beta[i,1] <- beta0[1]
	beta[i,2] <- beta0[2]
	omega[i] <- omega0
}

# Trace Plots
par(mfrow=c(2,2))
plot(beta[,1])
plot(beta[,2])
plot(omega)

# Examine Fit
beta1 <- mean(beta[1001:5000,1])
beta2 <- mean(beta[1001:5000,2])

plot(gdpgrowth$DEF60,gdpgrowth$GR6096,main='Fit of Bayesian Linear Model',xlab='Defense Spending as a Fraction of GDP',ylab='GDP Growth Rate')
abline(beta1,beta2)

# Heavy Tailed Error Model
X <- cbind(1,gdpgrowth$DEF60)
Y <- gdpgrowth$GR6096

lambda0 <- diag(1,dim(X)[1])
K <- diag(.1,2)
m <- matrix(0,2,1)
d <- 1
eta <- 1
omega0 <- 1
h <- 1


n.iter <- 5000

beta <- matrix(0,n.iter,2)
omega <- rep(0,n.iter)
lambda <- NULL

for (i in 1:n.iter){
	mu.star <- t((t(Y)%*%lambda0%*%X)%*%solve(t(X)%*%lambda0%*%X + K))
	lambda.star <- t(X)%*%lambda0%*%X + K
	eta.star <- eta + t(Y)%*%lambda0%*%Y-t(Y)%*%lambda0%*%X%*%solve(t(X)%*%lambda0%*%X + K)%*%t(t(Y)%*%lambda0%*%X)
	beta0 <- mvrnorm(1,mu.star,solve(omega0*lambda.star))
	omega0 <- rgamma(1,(d+dim(X)[1])/2,eta.star/2)
	beta[i,1] <- beta0[1]
	beta[i,2] <- beta0[2]
	omega[i] <- omega0
	for (j in 1:dim(X)[1]){
		lambda0[j,j] <- rgamma(1,(h+1)/2,.5*(Y[j]-t(X[j,])%*%beta0)^2+h)
	}
	lambda[[i]] <-lambda0 
}

beta11 <- mean(beta[1001:5000,1])
beta22 <- mean(beta[1001:5000,2])
plot(gdpgrowth$DEF60,gdpgrowth$GR6096,main='Fit of Heavy Tailed Error Model',xlab='Defense Spending as a Fraction of GDP',ylab='GDP Growth Rate')
abline(beta11,beta22)
abline(beta1,beta2,col=2)
