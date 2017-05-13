#######################################
##### Weighting & Kernal Function #####
#######################################

#Gaussian Kernal

k2 <- function(x){
	exp(-(x^2)/2)
}

#Weighting Function

w <- function(x,x.star,h,K){
	(1/h)*K((x-x.star)/h)
}

######################################
########## Cross Validation ##########
######################################

# Cross Validation: one training set, one test set
# Inputs: test data, training data, bandwith, and Kernal

cv <- function(test.x,test.y,train.x,train.y,h,K){
	pred.err <- matrix(0,length(h),2)

	for (j in 1:length(h)){
		xtest  <- test.x-mean(train.x)
		ytest  <- test.y-mean(train.y)
		xtrain <- train.x-mean(train.x)
		ytrain <- train.y-mean(train.y) 

		est <- rep(0,length(xtest))
		for (i in 1:length(xtest)){
			temp <- w(xtrain,xtest[i],h[j],K)
			est[i] <- sum((temp/sum(temp))*ytrain)
		}
		pred.err[j,1] <- h[j]
		pred.err[j,2] <- sum((ytest-est)^2)
	}
	pred.err
}

# N-fold Cross Validation
# Inputs: data, sample sive, number of folds, bandwith, Kernal

cv.N <- function(n,N,X,Y,h,K){	
	pred.err <- matrix(0,length(h),N)
	
	draws <- NULL
	for (i in 1:N){
		draws <- c(draws,rep(i,floor(n/N)))
	}
	draws <- c(draws,sample(1:N,n-N*floor(n/N)))
	samp  <- sample(draws,n)
	
	for (s in 1:N){
		train.y <- Y[-which(samp==s)]
		test.y  <- Y[which(samp==s)]
		train.x <- X[-which(samp==s)]
		test.x  <- X[which(samp==s)]
		
		for (j in 1:length(h)){	
			xtest  <- test.x-mean(train.x)
			ytest  <- test.y-mean(train.y)
			xtrain <- train.x-mean(train.x)
			ytrain <- train.y-mean(train.y) 

			est <- rep(0,length(xtest))
			for (i in 1:length(xtest)){
				temp <- w(xtrain,xtest[i],h[j],K)
				est[i] <- sum((temp/sum(temp))*ytrain)
			}
			pred.err[j,s] <- (1/length(test.x))*sum((ytest-est)^2)
		}
	}
	pred.err
}

#######################################
############ Generate Data ############
#######################################
# Functions

f1 <- function(x){
	x^2
}

f2 <- function(x){
	sin(4*pi*x)
}

#Generate data

x <- runif(500)

y1 <- f1(x)+rnorm(500,0,.2)
y2 <- f1(x)+rnorm(500,0,.5)
y3 <- f2(x)+rnorm(50,0,.2)
y4 <- f2(x)+rnorm(500,0,.5)

########################################
######### Run Cross Validation #########
########################################
h <- seq(0.001,0.5,by=.001)

# R cv function
samp <- sample(500,100)

train1y <- y1[-samp]
test1y  <- y1[samp]

train1x <- x[-samp]
test1x  <- x[samp]

ex1 <- cv(test1x,test1y,train1x,train1y,h,k2)

plot(ex1)
h[which.min(ex1[,2])]

# Run N-fold Cross Validation
n <- length(x)

plotcv.N <- function(n,N,X,Y,h,K){
ex <- cv.N(n,N,X,Y,h,K)
plot(h,ex[,1],type='l',ylim=c(min(ex),max(ex)),xlab='Bandwith',ylab='Prediction Error')
	for (i in 2:N){
		points(h,ex[,i],type='l',col=i)
	}
	h[which.min(rowSums(ex))]
}

par(mfrow=c(2,2))

h <- seq(0.001,0.5,by=.001)

plotcv.N(n,5,x,y1,h,k2)
# 0.036
plotcv.N(n,5,x,y2,h,k2)
# 0.095

h <- seq(0.001,0.1,by=.001)

plotcv.N(n,5,x,y3,h,k2)
# 0.015
plotcv.N(n,5,x,y4,h,k2)
# 0.025

########################################
############# Plots of Fit #############
########################################

xa  <- x - mean(x)
y1a <- y1 - mean(y1)
y2a <- y2 - mean(y2)
y3a <- y3 - mean(y3)
y4a <- y4 - mean(y4)

#Plot f1 - sigma=.2

plot(f1,main='Smooth - Low Noise',ylim=c(-1.2,1.75))
points(x,y1)

h <- c(0.036)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.01)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k2)
		x.est[i] <- sum((temp/sum(temp))*y1a)
	}

points(new.x+mean(x),x.est+mean(y1),type='l',col=j+1)
}

#Plot f1  - sigma=.5

plot(f1,main='Smooth - High Noise',ylim=c(-1.2,1.75))
points(x,y2)

h <- c(0.095)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.01)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k2)
		x.est[i] <- sum((temp/sum(temp))*y2a)
	}

points(new.x+mean(x),x.est+mean(y2),type='l',col=j+1)
}

#Plot f2 sigma=.2

plot(f2,main='Wiggly - Low Noise',ylim=c(-2.25,2.25))
points(x,y3)

h <- c(.015)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.01)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k2)
		x.est[i] <- sum((temp/sum(temp))*y3a)
	}

points(new.x+mean(x),x.est+mean(y3),type='l',col=j+1)
}

#Plot f2 sigma=.5

plot(f2,main='Wiggly - High Noise',ylim=c(-2.25,2.25))
points(x,y4)

h <- c(.025)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.01)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k2)
		x.est[i] <- sum((temp/sum(temp))*y4a)
	}

points(new.x+mean(x),x.est+mean(y4),type='l',col=j+1)
}

