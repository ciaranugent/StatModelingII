tests <- read.csv('https://raw.githubusercontent.com/jgscott/SDS383D/master/data/mathtest.csv')

ybar <- mean(tests[,2])

schools <- unique(tests[,1])

b <- matrix(0,length(schools),3)
for (i in 1:length(schools)){
	b[i,1] <- i
	b[i,2] <- mean(tests[tests[,1]==i,2])
	b[i,3] <- length(tests[tests[,1]==i,2])
}

# Set hyperparameters
a1 <- 1
b1 <- 1
a2 <- 1
b2 <- 1

mu0    <- ybar
sigma0 <- 10

n <- dim(tests)[1]
p <- length(schools)

# Set initial values
mu    <- ybar
theta <- b[,2]
sigma <- 1
tau   <- 1

# Gibbs Sampler

n.iter <- 10000

post.vars <- matrix(0,n.iter,103)

for (j in 1:n.iter){
	theta.vec <- (theta[tests[,1]])
	sigma <- 1/rgamma(1,(n+p)/2 + a1, sum((theta.vec-tests[,2])^2)/2 + sum((theta-mu)^2)/(2*tau) + b1)
	tau   <- 1/rgamma(1,p/2 + a2,sum((theta-mu)^2)/(2*sigma) + b2)
	mu    <- rnorm(1,(sum(theta)/(tau*sigma) + ybar/sigma0)/(p/(tau*sigma)+1/sigma0),sqrt(1/(p/(tau*sigma)+1/sigma0)))
	for (i in 1:p){
		theta[i] <- rnorm(1,(sum(tests[tests[,1]==i,2])/sigma + mu/(tau*sigma))/(b[i,3]/sigma +1/(tau*sigma)),sqrt(1/(b[i,3]/sigma +1/(tau*sigma))))
		post.vars[j,i] <- theta[i]
	} 
	post.vars[j,101] <- sigma
	post.vars[j,102] <- tau
	post.vars[j,103] <- mu
}

post.theta <- apply(post.vars[500:1000,1:100],2,mean)
sqrt(mean(post.vars[500:10000,101]))
#9.190537
sqrt(mean(post.vars[500:10000,102]))
#0.5499993
mean(post.vars[500:10000,103])
#48.10882

plot(b[,2],post.theta,ylab='Posterior Mean',xlab='Average by School')
abline(0,1)

par(mfrow=c(2,2))
plot(post.vars[500:10000,1],main='theta1',xlab='Iteration',ylab='theta1')
plot(post.vars[500:10000,101],main='sigma^2',xlab='Iteration',ylab='sigma')
plot(post.vars[500:10000,102],main='tau^2',xlab='Iteration',ylab='tau')
plot(post.vars[500:10000,103],main='mu',xlab='Iteration',ylab='mu')

hist(post.vars[500:10000,1],main='theta1',breaks=20,freq=FALSE,xlab='Value')
hist(post.vars[500:10000,101],main='sigma^2',breaks=20,freq=FALSE,xlab='Value')
hist(post.vars[500:10000,102],main='tau^2',breaks=20,freq=FALSE,xlab='Value')
hist(post.vars[500:10000,103],main='mu',breaks=20,freq=FALSE,xlab='Value')

k <- abs(ybar-post.theta)/ybar
plot(b[,3],k,xlab='Sample Size',ylab='K',main='Plot of Shrinkage by Sample Size')