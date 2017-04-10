library(lme4)
library(lattice)
library(MCMCpack)
library(MASS)
library(Matrix)

cheese  <- read.csv("https://raw.githubusercontent.com/jgscott/SDS383D/master/data/cheese.csv")

cheese0 <- cheese[cheese[,4]==0,]
cheese1 <- cheese[cheese[,4]==1,]

model.lm0 <- lm(log(cheese0$vol) ~ 1 + log(cheese0$price), data=cheese0)
model.lm1 <- lm(log(cheese1$vol) ~ 1 + log(cheese1$price), data=cheese1)

# Exploratory plots

xyplot( log(cheese$vol) ~ log(cheese$price) , group=as.factor(cheese$disp), 
       xlab='log(Price)', ylab='log(Volume)',main='Scatterplots of Data by Display',
       scales = list(x = list(log = 10, equispaced.log = FALSE)),
       type = c("p", "r"))

xyplot( log(cheese$vol) ~ log(cheese$price) | cheese$store, group=as.factor(cheese$disp), 
       strip=FALSE, xlab='log(Price)', ylab='log(Volume)',main='Scatterplots by Store',
       scales = list(x = list(log = 10, equispaced.log = FALSE)),
       type = c("p", "r"))

stores <- unique(cheese[,1])

stores.avg <- matrix(0,88,2)
for (i in 1:88){
	temp <- cheese[cheese[,1]==stores[i],]
	stores.avg[i,1] <- mean(log(temp[temp$disp==0,2]*temp[temp$disp==0,3]))
	stores.avg[i,2] <- mean(log(temp[temp$disp==1,2]*temp[temp$disp==1,3]))
}

boxplot(stores.avg,xlab='Display',ylab='Mean log(Sales)',main='Boxplot of Average Store Sales by Display',names=c('0','1'))

# Fitting the model using LMER

model.rml <- lmer(log(cheese$vol) ~ log(cheese$price) + cheese$disp + log(cheese$price):cheese$disp | store, data=cheese)

plot(model.rml)
qqnorm(residuals(model.rml))

coef(model.rml)

cheese28 <- cheese[cheese[,1]==stores[28],]
cheese88 <- cheese[cheese[,1]==stores[88],]

par(mfrow=c(1,2))
plot(log(cheese28$price),log(cheese28$vol),col=cheese28$disp+1,main='ATLANTA - WINN DIXIE',xlab='log(price)',ylab='log(volume)')
abline(8.049699,-0.05569669)
abline(8.049699-0.29894912,-0.05569669 + 0.22842926,col=2)

plot(log(cheese88$price),log(cheese88$vol),col=cheese88$disp+1,main='DALLAS/FT. WORTH - WINN DIXIE',xlab='log(price)',ylab='log(volume)')
abline(9.994213,-3.32491911)
abline(9.994213+1.56628926 ,-3.32491911-1.06414002,col=2)

# Heirarchical Bayesian Model via Gibbs Sampling

st <- length(stores)
n.iter <- 5000

X <- cbind(1,log(cheese[,2]),cheese[,4],log(cheese[,2])*cheese[,4])
Y <- log(cheese[,3])

# Set Hyper-Priors
Psi <- diag(rep(1,4))

lm(log(cheese$vol) ~ 1 + log(cheese$price) + cheese$disp + log(cheese$price):cheese$disp)$coefficients
#(Intercept)             log(cheese$price)                   cheese$disp log(cheese$price):cheese$disp 
#8.8095803                    -0.8898247                     0.7695531                    -0.3958158 
Theta.0 <- c(8.8095803, -0.8898247, 0.7695531, -0.3958158)

nu <- 4
p  <- 4

# Initialize values

SigmaInv <- solve(Psi)
Theta.i  <- Theta.0

beta<-NULL
beta.i <- NULL
for (i in 1:st){
	beta[[i]] <- matrix(0,n.iter,4)
	beta.i[[i]] <- rep(0,4)
}
sig.i <- 1
Sigma <- NULL


for(i in 1:n.iter){
betasum1 <- matrix(0,4,4)
totsum   <- 0
betasum2 <- matrix(0,4,1)
for (j in 1:st){
	temp.X <- X[cheese[,1]==stores[j],]
	temp.Y <- Y[cheese[,1]==stores[j]]
	xx <- matrix(0,4,4)
	yx <- matrix(0,4,1)
	for(k in 1:dim(temp.X)[1]){
		xx <- xx + temp.X[k,]%*%t(temp.X[k,])
		yx <- yx + temp.X[k,]*temp.Y[k]
	}	
	beta.i[[j]] <- mvrnorm(1,ginv(SigmaInv + sig.i*xx)%*%(SigmaInv%*%Theta.i + sig.i*yx), ginv(SigmaInv + sig.i*xx))
	beta[[j]][i,1] <- beta.i[[j]][1]
	beta[[j]][i,2] <- beta.i[[j]][2]
	beta[[j]][i,3] <- beta.i[[j]][3]
	beta[[j]][i,4] <- beta.i[[j]][4]
	# We will need these later
	betasum1 <- betasum1 + (beta.i[[j]] - Theta.i)%*%t(beta.i[[j]] - Theta.i)
	totsum   <- totsum + sum((temp.X%*%beta.i[[j]]-temp.Y)^2)
	betasum2 <- betasum2 + beta.i[[j]]
}

Sigma[[i]] <- riwish(nu+88,Psi+betasum1%*%t(betasum1))
SigmaInv <- ginv(Sigma[[i]])

sig.i   <- 1/rgamma(1,dim(X)[2]/2,totsum/2)
Theta.i <- mvrnorm(1,ginv(diag(rep(1,4))+88*SigmaInv)%*%(Theta.0 + SigmaInv%*%betasum2),ginv(diag(rep(1,4))+88*SigmaInv))
}

# Plots of Results for specific stores

cheese28 <- cheese[cheese[,1]==stores[28],]
cheese88 <- cheese[cheese[,1]==stores[88],]

par(mfrow=c(1,2))
plot(log(cheese28$price),log(cheese28$vol),col=cheese28$disp+1,main='ATLANTA - WINN DIXIE',xlab='log(price)',ylab='log(volume)')
abline(colMeans(beta[[28]])[1],colMeans(beta[[28]])[2])
abline(colMeans(beta[[28]])[1]+colMeans(beta[[28]])[3],colMeans(beta[[28]])[2]+colMeans(beta[[28]])[4],col=2)

plot(log(cheese88$price),log(cheese88$vol),col=cheese88$disp+1,main='DALLAS/FT. WORTH - WINN DIXIE',xlab='log(price)',ylab='log(volume)')
abline(colMeans(beta[[88]])[1],colMeans(beta[[88]])[2])
abline(colMeans(beta[[88]])[1]+colMeans(beta[[88]])[3],colMeans(beta[[88]])[2]+colMeans(beta[[88]])[4],col=2)

