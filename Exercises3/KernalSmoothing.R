#######################################
############ Generate Data ############
#######################################

# Functions

f1 <- function(x){
	x^2
}

f2 <- function(x){
	sin(2*pi*x)
}

#Generate data

x <- runif(50)

y1 <- f1(x)+rnorm(50,0,.2)
y2 <- f2(x)+rnorm(50,0,.2)

# Center data

y1a <- y1 -mean(y1)
y2a <- y2 - mean(y2)
xa  <- x - mean(x)

#######################################
###### Build Weighting Functions ######
#######################################

#Kernal Functions

#Uniform Kernal

k1 <- function(x){
	as.numeric(ifelse(abs(x)<=1,1,0))
}

#Gaussian Kernal

k2 <- function(x){
	exp(-(x^2)/2)
}

#Weighting Function

w <- function(x,x.star,h,K){
	(1/h)*K((x-x.star)/h)
}

#######################################
################ Plots ################
#######################################

#Plot f1 - Gaussian Kernal

plot(f1,main='Gaussian Kernal')
points(x,y1)

h <- c(.05,.1,.25,.5)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.01)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k2)
		x.est[i] <- sum((temp/sum(temp))*y1a)
	}

points(new.x+mean(x),x.est+mean(y1),type='l',col=j+1)
}
legend(.05,.95,legend=c('f(x)','h=.05','h=.1','h=.25','h=.5'),lty=1,col=1:5)

#Plot f1 - Uniform Kernal

plot(f1,main='Uniform Kernal')
points(x,y1)

h <- c(.1,.25,.5,.75)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.001)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k1)
		x.est[i] <- sum((temp/sum(temp))*y1a)
	}

points(new.x+mean(x),x.est+mean(y1),type='l',col=j+1)
}
legend(.05,.95,legend=c('f(x)','h=.1','h=.25','h=.5','h=.75'),lty=1,col=1:5)

#Plot f2 Gaussian Kernal

plot(f2,main='Gaussian Kernal',ylim=c(-1.1,1.1))
points(x,y2)

h <- c(.05,.1,.25,.5)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.01)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k2)
		x.est[i] <- sum((temp/sum(temp))*y2a)
	}

points(new.x+mean(x),x.est+mean(y2),type='l',col=j+1)
}
legend(.75,1,legend=c('f(x)','h=.05','h=.1','h=.25','h=.5'),lty=1,col=1:5)

#Plot f2 Uniform Kernal

plot(f2,main='Uniform Kernal',ylim=c(-1.1,1.1))
points(x,y2)

h <- c(.05,.1,.25,.5)

for (j in 1:length(h)){
new.x <- seq(0,1,by=.001)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k1)
		x.est[i] <- sum((temp/sum(temp))*y2a)
	}

points(new.x+mean(x),x.est+mean(y2),type='l',col=j+1)
}
legend(.75,1,legend=c('f(x)','h=.05','h=.1','h=.25','h=.5'),lty=1,col=1:5)

######################################
########## Cross Validation ##########
######################################
h <- seq(0:.5,by=.01)



for (j in 1:length(h)){
new.x <- seq(0,1,by=.001)-mean(x)

x.est <- rep(0,length(new.x))
	for (i in 1:length(new.x)){
		temp <- w(xa,new.x[i],h[j],k1)
		x.est[i] <- sum((temp/sum(temp))*y2a)
	}

points(new.x+mean(x),x.est+mean(y2),type='l',col=j+1)
}
