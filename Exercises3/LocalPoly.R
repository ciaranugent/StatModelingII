utilities <- read.table(file="/users/ciaranugent/Documents/utilities.txt",sep=',',header=TRUE)

temp <- utilities[,1]
cost <- utilities[,2]/utilities[,3]

plot(temp,cost)

# Gaussian Kernal

k2 <- function(x){
	exp(-(x^2)/2)
}

# Weighting Function

w <- function(x,x.star,h,K){
	(1/h)*K((x-x.star)/h)
}

# Local polynomial regression function

locpoly <- function(X,Y,new,h,D,K){
	fit.new <- NULL
	for (j in 1:length(new)){
		W <- (1/h)*diag(w(X,new[j],h,K))
		R <- NULL
		for (i in 0:D){
			R <- cbind(R,(X-new[j])^(i))
		}
		fit.new <- c(fit.new,(solve(t(R)%*%W%*%R)%*%t(R)%*%W%*%Y)[1])
	}
	fit.new
}

h.hat <- function(X,Y,new,h,D,K){
	fit.new <- NULL
	for (j in 1:length(new)){
		W <- (1/h)*diag(w(X,new[j],h,K))
		R <- NULL
		for (i in 0:D){
			R <- cbind(R,(X-new[j])^(i))
		}
		fit.new <- c(fit.new,(solve(t(R)%*%W%*%R)%*%t(R)%*%W%*%Y)[1])
	}
	fit.new	
}

# Leave one out cross-validation function

LOOCV <- function(X,Y,prop,D,K){
	err <- NULL
	for (q in 1:length(prop)){	
		temp.w <- 0
		for (w in 1:length(X)){
			temp.w <- (locpoly(X[-w],Y[-w],X[w],prop[q],D,K) - Y[w])^2 + temp.w
		}
		err <- c(err,temp.w)
	}
	err
}

prop.h <- seq(1,10,.1)

err.h <- LOOCV(temp,cost,prop.h,0,k2)
prop.h[which.min(err.h)]
# 3.9

err.h <- LOOCV(temp,cost,prop.h,1,k2)
prop.h[which.min(err.h)]
# 6.9

err.h <- LOOCV(temp,cost,prop.h,2,k2)
prop.h[which.min(err.h)]


new.x <- seq(10,80,.25)

# Plots
plot(temp,cost,main='Plot of Local Polynomial Fits')
new.y <- locpoly(temp,cost,new.x,3.9,0,k2)
points(new.x,new.y,type='l')
new.y <- locpoly(temp,cost,new.x,6.9,1,k2)
points(new.x,new.y,type='l',col=2)
new.y <- locpoly(temp,cost,new.x,90.9,2,k2)
points(new.x,new.y,type='l',col=3)
legend(60,7,legend=c('Intercept','Linear','Quadratic'),lty=1,col=1:3)

y.fit <- locpoly(temp,cost,temp,6.9,1,k2)
RSS <- sum((cost-y.fit)^2)/(length(temp)-1)

plot(temp,cost-y.fit,ylab='Residuals',xlab='Temperature',main='Plot of Residuals for Local Linear Regression')


# log transformation

err.h <- LOOCV(temp,log(cost),prop.h,1,k2)
prop.h[which.min(err.h)]
#5.4

y.fit <- locpoly(temp,log(cost),temp,6.9,1,k2)
RSS <- sum((log(cost)-y.fit)^2)/(length(temp)-1)
plot(temp,log(cost)-y.fit,ylab='Residuals',xlab='Temperature',main='Plot of Residuals - Log Transformed Cost')

plot(temp,log(cost),main='Plot of Local Polynomial Fit')
new.y <- locpoly(temp,log(cost),new.x,5.4,1,k2)
points(new.x,new.y,type='l')
points(new.x,new.y+1.96*RSS,type='l',col=2,lty=2)
points(new.x,new.y-1.96*RSS,type='l',col=2,lty=2)



