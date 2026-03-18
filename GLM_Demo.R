setwd('C:/Users/chaslebacher/Desktop/R_code/Europa_Bingham')

source("Directional.R")
source("GLM_3D.R")

# Test Sample2_vMF() ----
n0 <- 100
Z <- matrix(0,3*n0,3)
Z[1:(2*n0),1] <- 3
Z[1:(2*n0),2] <- -1
Z[(2*n0+1):(3*n0),2] <- 5
col <- c(rep('blue',200),rep('red',100))
Y <- Sample2_vMF(Z)
plot(Y[,1],Y[,2],col=col,pch=19,xlim=c(-1,1),ylim=c(-1,1))


# Illustrate GLM_vMF ----

## Multiple linear regression ----

# Generate covariates:
r0 <- 4
m0 <- 20
n <- m0*(2*r0 + 1)^2
X <- matrix(0,n,2)
X[,1] <- rep(rep((-r0):r0,each=2*r0+1)/r0,each=m0)
X[,2] <- rep(rep((-r0):r0,times=2*r0+1)/r0,each=m0)
plot(X[,1],X[,2])
# True parameter matrix:
Theta0 <- matrix(c(1,0,1,0,0,2),2,3)
Theta0
Z0 <- cbind(1,X) %*% t(Theta0)
Z0[sample(n,10),]
# Generate responses:
Y <- Sample2_vMF(Z0)
str(Y)

k1 <- cos(2*pi*seq(0,1,0.01))
k2 <- sin(2*pi*seq(0,1,0.01))
plot(X[,1],X[,2],pch=16,col='gray',
	 xlim=c(-r0-1,r0+1)/r0,ylim=c(-r0-1,r0+1)/r0)
for (x in (-r0):r0){
	for (y in (-r0):r0){
		lines((x+k1/2)/r0,(y+k2/2)/r0,col='gray')
	}
}
for (x in (-r0):r0){
	for (y in (-r0):r0){
		lines((x+k1/2)/r0,(y+k2/2)/r0,col='gray')
		z <- Theta0 %*% c(1,x/r0,y/r0)
		mz <- Mu_vMF(z)
		lines(c(x,x+mz[1]/2)/r0,c(y,y+mz[2]/2)/r0,
			  col='green',lwd=4)
	}
}
for (i in 1:n){
	lines(X[i,1] + c(0,Y[i,1]/2/r0),
		  X[i,2] + c(0,Y[i,2]/2/r0),col='gray')
}
for (x in (-r0):r0){
	for (y in (-r0):r0){
		mz <- colMeans(Y[X[,1]==x/r0 & X[,2]==y/r0,])
		lines(c(x,x+mz[1]/2)/r0,c(y,y+mz[2]/2)/r0,
			  col='black',lwd=1)
	}
}


### Fit the GLM ----
Theta <- GLM_vMF(Y,X,W=rep(1,n))
Theta
Theta0
Theta - Theta0

plot(X[,1],X[,2],pch=16,col='gray',
	 xlim=c(-r0-1,r0+1)/r0,ylim=c(-r0-1,r0+1)/r0)
for (x in (-r0):r0){
	for (y in (-r0):r0){
		lines((x+k1/2)/r0,(y+k2/2)/r0,col='gray')
	}
}
for (x in (-r0):r0){
	for (y in (-r0):r0){
		lines((x+k1/2)/r0,(y+k2/2)/r0,col='gray')
		z <- Theta0 %*% c(1,x/r0,y/r0)
		mz <- Mu_vMF(z)
		lines(c(x,x+mz[1]/2)/r0,c(y,y+mz[2]/2)/r0,
			  col='green',lwd=4)
	}
}
for (x in (-r0):r0){
	for (y in (-r0):r0){
		lines((x+k1/2)/r0,(y+k2/2)/r0,col='gray')
		mz <- colMeans(Y[X[,1]==x/r0 & X[,2]==y/r0,])
		lines(c(x,x+mz[1]/2)/r0,c(y,y+mz[2]/2)/r0,
			  col='black',lwd=1)
	}
}
for (x in (-r0):r0){
	for (y in (-r0):r0){
		z <- Theta %*% c(1,x/r0,y/r0)
		mz <- Mu_vMF(z)
		lines(c(x,x+mz[1]/2)/r0,c(y,y+mz[2]/2)/r0,
			  col='blue',lwd=2)
	}
}

### Root mean squared errors of estimated means ----
error_glm <- 0
error_emp <- 0
for (x in ((-r0):r0)/r0){
	for (y in ((-r0):r0)/r0){
		z0 <- Theta0 %*% c(1,x,y)
		mz0 <- Mu_vMF(z0)
		z_glm <- Theta %*% c(1,x,y)
		mz_glm <- Mu_vMF(z_glm)
		mz_emp <- colMeans(Y[X[,1]==x & X[,2]==y,])
		error_glm <- error_glm +
			sqrt(sum((mz_glm - mz0)^2))
		error_emp <- error_emp +
			sqrt(sum((mz_emp - mz0)^2))
	}
}
error_glm <- error_glm/(2*r0 + 1)^2
error_emp <- error_emp/(2*r0 + 1)^2
c('glm'=error_glm,'emp'=error_emp)


## Multiple quadratic regression ----

# Generate covariates:
r0 <- 10
n <- (2*r0 + 1)^2
X <- matrix(1,n,5)
colnames(X) <- c('X1','X2','X1^2','X1*X2','X2^2')
X[,1] <- rep((-r0):r0,each=2*r0+1)/r0
X[,2] <- rep((-r0):r0,times=2*r0+1)/r0
X[,3] <- X[,1]^2
X[,4] <- X[,1]*X[,2]
X[,5] <- X[,2]^2
plot(X[,1],X[,2])
# True parameter matrix:
Theta0 <- t(matrix(c(1,2,0,-2,0,0,-1,0,-1,0,1,-1),6,2))
colnames(Theta0) <- c('const',colnames(X))
Theta0
Z0 <- cbind(1,X) %*% t(Theta0)
Z0[sample(n,10),]

# Plot parameters:
plot(X[,1],X[,2],pch=16,col='gray',
	 xlim=c(-r0-1,r0+1)/r0,ylim=c(-r0-1,r0+1)/r0)
magnfac <- 0.04
for (i in 1:n){
	lines(X[i,1] + magnfac*c(0,Z0[i,1]),
		  X[i,2] + magnfac*c(0,Z0[i,2]),
		  col='forestgreen')
}

# Generate responses:
Y <- Sample2_vMF(Z0)
str(Y)

# Plot means and data:
k1 <- cos(2*pi*seq(0,1,0.01))
k2 <- sin(2*pi*seq(0,1,0.01))
plot(X[,1],X[,2],pch='.',col='gray',
	 xlim=c(-r0-1,r0+1)/r0,ylim=c(-r0-1,r0+1)/r0)
for (i in 1:n){
	lines(X[i,1] + k1/2/r0,X[i,2] + k2/2/r0,col='gray')
}
M <- matrix(0,n,2)
for (i in 1:n){
	z <- Theta0 %*% c(1,X[i,])
	mz <- Mu_vMF(z)
	M[i,] <- mz
	lines(X[i,1] + c(0,mz[1]/2/r0),
		  X[i,2] + c(0,mz[2]/2/r0),col='green',lwd=3)
	lines(X[i,1] + c(0,Y[i,1]/2)/r0,
		  X[i,2] + c(0,Y[i,2]/2)/r0,
		  col='black',lwd=1)
}

### Fit the GLM ----
Theta <- GLM_vMF(Y,X,W=rep(1,n))
Theta
Theta0
Theta - Theta0

Mhat <- matrix(0,n,2)
for (i in 1:n){
	z <- Theta %*% c(1,X[i,])
	mz <- Mu_vMF(z)
	Mhat[i,] <- mz
	lines(X[i,1] + c(0,mz[1]/2/r0),
		  X[i,2] + c(0,mz[2]/2/r0),col='blue',lwd=1)
}

# Plot true and estimated parameters:
plot(X[,1],X[,2],pch=16,col='gray',
	 xlim=c(-r0-1,r0+1)/r0,ylim=c(-r0-1,r0+1)/r0)
magnfac <- 0.04
Z <- cbind(1,X) %*% t(Theta)
for (i in 1:n){
	lines(X[i,1] + magnfac*c(0,Z0[i,1]),
		  X[i,2] + magnfac*c(0,Z0[i,2]),
		  col='green',lwd=2)
	lines(X[i,1] + magnfac*c(0,Z[i,1]),
		  X[i,2] + magnfac*c(0,Z[i,2]),
		  col='blue')
}

### Mean errors of estimated means ----
error_glm <- mean(sqrt(rowSums((Mhat - M)^2)))
error_emp <- mean(sqrt(rowSums((Y    - M)^2)))
c('glm'=error_glm,'emp'=error_emp)


### Fit the wrong GLM ----
X_red <- X[,1:2]
Theta_red <- GLM_vMF(Y,X_red,W=rep(1,n))
Theta_red

Mhat <- matrix(0,n,2)
for (i in 1:n){
	z <- Theta_red %*% c(1,X_red[i,])
	mz <- Mu_vMF(z)
	Mhat[i,] <- mz
	# lines(X[i,1] + c(0,mz[1]/2/r0),
	# 	  X[i,2] + c(0,mz[2]/2/r0),col='blue',lwd=1)
}

# Plot true and estimated parameters:
plot(X[,1],X[,2],pch=16,col='gray',
	 xlim=c(-r0-1,r0+1)/r0,ylim=c(-r0-1,r0+1)/r0)
magnfac <- 0.05
Z_red <- cbind(1,X_red) %*% t(Theta_red)
for (i in 1:n){
	lines(X[i,1] + magnfac*c(0,Z0[i,1]),
		  X[i,2] + magnfac*c(0,Z0[i,2]),
		  col='green',lwd=2)
	lines(X[i,1] + magnfac*c(0,Z_red[i,1]),
		  X[i,2] + magnfac*c(0,Z_red[i,2]),
		  col='blue')
}


## Test of target function and its derivatives ----

X1 <- cbind(1,X[,1:2])
Theta <- t(qr.solve(X1,Y))
Delta <- matrix(rnorm(2*3,0,0.5),2,3)
res0 <- TestGLM(Y,X1,Theta)
t0 <- 0.0001
res1 <- TestGLM(Y,X1,Theta + t0*Delta)
c(res0$FTheta, res1$FTheta)
c((res1$FTheta - res0$FTheta)/t0,sum(res0$Grad * Delta))
c(2*(res1$FTheta - res0$FTheta - t0*sum(res0$Grad * Delta))/t0^2,
  t(as.vector(Delta)) %*% res0$Hessian %*% as.vector(Delta))

(res1$Grad - res0$Grad)/t0 -
	matrix(res0$Hessian %*% as.vector(Delta),2,3)

Delta <- matrix(qr.solve(res0$Hessian,as.vector(res0$Grad)),2,3)
Delta

t0 <- 0.0001
res1 <- TestGLM(Y,X1,Theta - t0*Delta)
c(res0$FTheta, res1$FTheta)
c((res1$FTheta - res0$FTheta)/t0,-sum(res0$Grad * Delta))
c(2*(res1$FTheta - res0$FTheta + t0*sum(res0$Grad * Delta))/t0^2,
  t(as.vector(Delta)) %*% res0$Hessian %*% as.vector(Delta))

(res1$Grad - res0$Grad)/t0 +
	matrix(res0$Hessian %*% as.vector(Delta),2,3)


## Local polynomial models ----

# Generate a grid of points:
r0 <- 10
nn <- (2*r0 + 1)^2
xx <- matrix(0,nn,2)
xx[,1] <- rep((-r0):r0,each=2*r0+1)/r0
xx[,2] <- rep((-r0):r0,times=2*r0+1)/r0
plot(xx[,1],xx[,2])

# Generate covariates:
set.seed(3032024)
n <- 4000
X <- round(matrix(2*runif(2*n)-1,n,2),2)

# True parameter matrices:
a1 <- 2
a2 <- 1
a3 <- 3
Z0 <- exp(-a1*rowSums(X^2))*cbind(a2,a3*X[,1])
Z0xx <- exp(-a1*rowSums(xx^2))*cbind(a2,a3*xx[,1])
# True mean matrices:
M0 <- Z0
for (i in 1:n){
	M0[i,] <- Mu_vMF(M0[i,])
}
M0xx <- Z0xx
for (i in 1:nn){
	M0xx[i,] <- Mu_vMF(M0xx[i,])
}

# Visualize model:
VisualiseXYV(xx,Y=Z0xx,magfac=0.18)
VisualiseXYV(xx,Y=M0xx,magfac=0.18)

# Generate responses:
Y <- Sample2_vMF(Z0)
str(Y)

# Visualise raw data:
VisualiseXYV(X,Y=Y,magfac=0.05)
VisualiseXYV_local(X,Y=Y,N=100,magfac=0.05) -> x0
VisualiseXYV(X,Y=M0,magfac=0.15)

x0 <- c(-0.95,0.95)
x0 <- c(0.5,0)
VisualiseXYV_local(X,Y=Y,N=400,x0=x0,magfac=0.05)


# Fit the nonparametric models:
LPFxx <- LocalPolynGLM_vMF(Y,X,xx,N=400,regfac=0,
							   showsteps=TRUE)
str(LPFxx)
saveRDS(LPFxx,file='GLM_N400.rds')

# Visualise true and estimated parameters z:
Zhxx <- LPFxx$Zcxx
Zhxx <- LPFxx$Zlxx
Zhxx <- LPFxx$Zqxx
magfac <- 0.18
VisualiseXYV(xx,Z0xx,colY='green',lwdY=2,magfac=magfac)
for (j in 1:nn){
	lines(xx[j,1] + magfac*c(0,Zhxx[j,1]),
		  xx[j,2] + magfac*c(0,Zhxx[j,2]),
		  col='blue',lwd=1)
}

# Visualise true and estimated parameters mu(z):
Mhxx <- LPFxx$Mcxx
Mhxx <- LPFxx$Mlxx
Mhxx <- LPFxx$Mqxx
magfac <- 0.18
VisualiseXYV(xx,M0xx,colY='green',lwdY=3,magfac=magfac)
for (j in 1:nn){
	lines(xx[j,1] + magfac*c(0,M0xx[j,1]),
		  xx[j,2] + magfac*c(0,M0xx[j,2]),
		  col='forestgreen',lwd=1)
	lines(xx[j,1] + magfac*c(0,Mhxx[j,1]),
		  xx[j,2] + magfac*c(0,Mhxx[j,2]),
		  col='black',lwd=2)
}

# Error measures:
rmse <- rep(Inf,3)
names(rmse) <- c('const','lin','quadr')
rmse[1] <- sqrt(mean(rowSums((LPFxx$Mcxx - M0xx)^2)))
rmse[2] <- sqrt(mean(rowSums((LPFxx$Mlxx - M0xx)^2)))
rmse[3] <- sqrt(mean(rowSums((LPFxx$Mqxx - M0xx)^2)))
rmse
round(rmse,4)


## Axial densities ----

xy <- rbind(c(-0.5,-0.5),c(1.5,-0.5),c(-0.5,1.5),c(1.5,1.5))
Z <- rbind(c(0,1),c(1,1),c(1,-1),c(0,0))
plot(xy[,1],xy[,2],xlim=c(-1.5,2.5),ylim=c(-1.5,2.5),xlab='x[1]',ylab='x[2]')
PlotAxialDensity(xy=xy,Z=Z,magfac=0.5)
