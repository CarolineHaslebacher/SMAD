# Fitting a GLM ----

GLM_vMF <- function(Y,X,add.constant=TRUE,
					W=rep(1,dim(Y)[1]),
					Theta0=NULL, # potential starting parameter
					show.steps=TRUE,
					prec=10^(-11),regfac=0)
	# Multiple generalized linear regression for directional data
{
	# Initializing:
	n <- dim(Y)[1]
	d <- dim(Y)[2]
	if (add.constant){
		X1 <- cbind(1,X)
		if (is.null(colnames(X))){
			colnames(X1) <- c('const',paste('X',1:dim(X)[2],sep=''))
		}else{
			colnames(X1) <- c('const',colnames(X))
		}
	}else{
		X1 <- X
	}
	r <- dim(X1)[2]
	# Starting parameter:
	Theta <- d*t(qr.solve(X1,Y))
	# Target function:
	ZX1 <- X1 %*% t(Theta)
	S <- sqrt(rowSums(ZX1^2))
	FTheta <- sum(W*lG0(S,d,prec)) - sum(W*Y*ZX1)
	# Comparison with proposed starting parameter:
	if (!is.null(Theta0)){
		ZX10 <- X1 %*% t(Theta0)
		S0 <- sqrt(rowSums(ZX10^2))
		FTheta0 <- sum(W*lG0(S0,d,prec)) - sum(W*Y*ZX10)
		if (FTheta0 < FTheta){
			Theta <- Theta0
			ZX1 <- ZX10
		}
	}
	# Start of iteration:
	Grad <- matrix(0,d,r)
	Hessian <- matrix(0,d*r,d*r)
	lastimprovement <- Inf
	iter1 <- 0
	while (lastimprovement > prec & iter1 < 500){
		iter1 <- iter1 + 1
		# Compute gradient and Hessian:
		Grad[] <- 0
		Hessian[] <- 0
		for (i in 1:n){
			tmp <- Moments_vMF(ZX1[i,],d,prec)
			Grad <- Grad + W[i]*
				(tmp$muz - Y[i,]) %*% t(X1[i,])
			Hessian <- Hessian + W[i]*
				kronecker(X1[i,]%*%t(X1[i,]),tmp$Sigmaz)
		}
		# Regularize Hessian:
		Hessian <- Hessian + regfac*mean(diag(Hessian))*diag(d*r)
		# Compute new candidate:
		Thetanew <- Theta - matrix(qr.solve(Hessian,as.vector(Grad)),d,r)
		ZX1 <- X1 %*% t(Thetanew)
		S <- sqrt(rowSums(ZX1^2))
		FThetanew <- sum(W*lG0(S,d,prec)) - sum(W*Y*ZX1)
		iter0 <- 0
		while (FThetanew >= FTheta & iter0 < 50){
			iter0 <- iter0+1
			Thetanew <- (Theta + Thetanew)/2
			ZX1 <- X1 %*% t(Thetanew)
			S <- sqrt(rowSums(ZX1^2))
			FThetanew <- sum(W*lG0(S,d,prec)) - sum(W*Y*ZX1)
		}
		if (iter0 == 50){
			lastimprovement <- 0 
		}else{
			lastimprovement <- FTheta - FThetanew
			Theta <- Thetanew
			FTheta <- FThetanew
			if (show.steps){
				print(round(c(iter1,FThetanew),5))
			}
		}
	}
	return(Theta)
}



# To test gradient and Hessian:

TestGLM <- function(Y,X1,Theta,prec=10^(-11))
{
	# Initializing:
	n <- dim(Y)[1]
	d <- dim(Y)[2]
	r <- dim(X1)[2]
	ZX <- X1 %*% t(Theta)
	S <- sqrt(rowSums(ZX^2))
	FTheta <- sum(lG0(S,d,prec)) - sum(Y*ZX)
	Grad <- matrix(0,d,r)
	Hessian <- matrix(0,d*r,d*r)
	for (i in 1:n){
		tmp <- Moments_vMF(ZX[i,],d,prec)
		Grad <- Grad +
			(tmp$muz - Y[i,]) %*% t(X1[i,])
		Hessian <- Hessian +
			kronecker(X1[i,]%*%t(X1[i,]),tmp$Sigmaz)
	}
	return(list(FTheta=FTheta,Grad=Grad,Hessian=Hessian))
}


# Local polynomial modelling ----

LocalPolynGLM_vMF <- function(Y,X,xx=X,
								 xsf=NULL,
								 bandwidth=1,N=NULL,regfac=0,
								 showsteps=FALSE, customweights=NULL)
	# Locally constant, linear and quadratic generalized linear
	# models for directional data with Euclidean covariate vector.
	# Output are the following tables:
	# -  xx
	# -  Zcxx, Zlxx, Zqxx:
	#    The i-th row of these matrices contains the estimated
	#    parameter for location xx[i,] via locally constant,
	#    linear and quadratic modelling, respectively.
	# -  Mcxx, Mlxx, Mqxx:
	#    The rows of these tables contain the corresponding
	#    estimated means.
{
	n <- dim(Y)[1]
	d <- dim(Y)[2]
	X <- as.matrix(X)
	q <- dim(X)[2]
	if (is.vector(xx)){
		xx <- matrix(xx,1,d)
	}
	nn <- dim(xx)[1]
	# Matrices of estimated vMF-parameters and means
	# print(d)
	# print(nn)
	# at locations in xx:
	Zcxx <- matrix(0,nn,d)
	Mcxx <- matrix(0,nn,d)
	Zlxx <- matrix(0,nn,d)
	Mlxx <- matrix(0,nn,d)
	Zqxx <- matrix(0,nn,d)
	Mqxx <- matrix(0,nn,d)
	
	# Scale factors for the columns of X and xx:
	if (is.null(xsf)){
		xsf <- apply(X,2,sd)
	}
	# Design matrix for GLM_vMF:
	DD <- matrix(1,n,(q+1)*(q+2)/2)
	# Parameter matrix (also as starting value):
	for (j in 1:nn){
		Xc <- t((t(X) - xx[j,])/xsf)
		DD[,2:(q+1)] <- Xc
		k <- q+1
		for (ell in 1:q){
			for (m in ell:q){
				k <- k+1
				DD[,k] <- Xc[,ell]*Xc[,m]
			}
		}
		if (is.null(N)){
			W <- exp(-rowSums(Xc^2)/2/bandwidth^2)
		}else{
			if (n > N){
				W <- SetWeights(rowSums(Xc^2),N) # CH: rowSums(Xc^2) is the distance 
			}else{
				W <- rep(1,n)
			}
		}
		# add customweights by multiplication: 
		if (!is.null(customweights)){
			W <- W*customweights
		}		
		
		# Locally constant fitting:
		Mcxx[j,] <- colSums(Y*W)/sum(W)
		Zcxx[j,] <- MLE0_vMF(Mcxx[j,])
		
		# Locally linear fitting:
		Theta <- GLM_vMF(Y,DD[,1:(q+1)],add.constant=FALSE,W=W,
						 Theta0=cbind(Zcxx[j,],
						 			 matrix(0,d,q)),
						 show.steps=FALSE,regfac=regfac)
		Zlxx[j,] <- Theta[,1]
		Mlxx[j,] <- Moments_vMF(Theta[,1])$muz
		
		# Locally quadratic fitting:
		Theta <- GLM_vMF(Y,DD,add.constant=FALSE,W=W,
					Theta0=cbind(Theta,matrix(0,d,q*(q+1)/2)),
					show.steps=FALSE,regfac=regfac)
		Zqxx[j,] <- Theta[,1]
		Mqxx[j,] <- Moments_vMF(Theta[,1])$muz
		if (showsteps){
			print(round(c(j,'x'=xx[j,],
						  'Zc'=Zcxx[j,],
						  'Zl'=Zlxx[j,],
						  'Zq'=Zqxx[j,]),3))
		}
	}
	return(list('xx'=xx,
				'Zcxx'= Zcxx,'Mcxx'=Mcxx,
				'Zlxx'= Zlxx,'Mlxx'=Mlxx,
				'Zqxx'= Zqxx,'Mqxx'=Mqxx))
}

SetWeights <- function(S,N,prec=10^(-11))
	# For a nonzero vector S of nonnegative numbers
	# and a target value N in (0,length(S)), this procedure
	# finds a value gamma > 0 such that the vector
	#    W = exp(-gamma*S)
	# satisfies sum(W) = N, and it returns W.
{
	gamma <- 0
	W <- rep(1,length(S))
	Hgamma <- length(S)
	gamma_new <- gamma + (Hgamma - N)/sum(S*W)
	while (gamma_new > (1 + prec)*gamma){
		gamma <- gamma_new
		W <- exp(-gamma*S)
		Hgamma <- sum(W)
		gamma_new <- gamma + (Hgamma - N)/sum(S*W)
	}
	return(W)
}

SunWeights <- function(S,N,prec=10^(-11))
	# For a nonzero vector S of nonnegative numbers
	# and a target value N in (0,length(S)), this procedure
	# finds a value gamma > 0 such that the vector
	#    W = exp(-gamma*S)
	# satisfies sum(W) = N, and it returns W.
{
	gamma <- 0
	W <- rep(1,length(S))
	Hgamma <- length(S)
	gamma_new <- gamma + (Hgamma - N)/sum(S*W)
	while (gamma_new > (1 + prec)*gamma){
		gamma <- gamma_new
		W <- exp(-gamma*S)
		Hgamma <- sum(W)
		gamma_new <- gamma + (Hgamma - N)/sum(S*W)
	}
	return(W)
}


# Simulate a response vector from vMF distributions ----

Sample2_vMF <- function(Z)
	# Z is a data matrix of n d-dimensional row vectors,
	# and the procedure simulates a response data matrix Y
	# such that Y[i,] follows the vMF distribution with
	# parameter vector Z[i,]
{
	n <- dim(Z)[1]
	d <- dim(Z)[2]
	S <- sqrt(rowSums(Z^2))
	Y <- matrix(0,n,d)
	notyet <- rep(TRUE,n)
	while (any(notyet)){
		for (i in (1:n)[notyet]){
			y <- rnorm(d)
			y <- y/sqrt(sum(y^2))
			if (runif(1) <= exp(sum(y*Z[i,]) - S[i])){
				Y[i,] <- y
				notyet[i] <- FALSE
			}
		}
	}
	return(Y)
}

# Visualise raw directional or axial data ----

VisualiseXYV <- function(X,Y=NULL,V=NULL,magfac=0.5,
						 pchX=16,colX='gray',
						 lwdY=1,colY='black',
						 lwdV=1,colV='blue',
						 xlim=NULL,ylim=NULL)
	# Input are data matrices
	# -  X (n x 2) : containing locations X[i,],
	# -  Y (n x 2) : containing direction vectors Y[i,] at X[i,]
	#                (not necessarily unit vectors),
	# -  V (n x 2) : containing axis vectors V[i,] (and
	#                implicitly -V[i,]) at X[i,]
	#                (not necessarily unit vectors).
	# The procedure plots the locations XX[i,] and the
	# corresponding directions Y[i,] and axes V[i,].
{
	n <- dim(X)[1]
	if (is.null(xlim)){
		xlim <- range(X[,1])
	}
	if (is.null(ylim)){
		ylim <- range(X[,2])
	}
	par(cex=1.2,mai=c(0.7,0.75,0.1,0.1),mgp=c(1.7,0.5,0))
	plot(X[,1],X[,2],pch=pchX,col=colX,
		 xlab=expression(italic(x[1])),xlim=xlim,
		 ylab=expression(italic(x[2])),ylim=ylim,
		 main='')
	if (!is.null(Y)){
		for (i in (1:n)){
			lines(X[i,1] + magfac*Y[i,1]*c(0,1),
				  X[i,2] + magfac*Y[i,2]*c(0,1),
				  lwd=lwdY,col=colY)
		}
	}
	if (!is.null(V)){
		for (i in (1:n)){
			lines(X[i,1] + magfac*V[i,1]*c(-1,1),
				  X[i,2] + magfac*V[i,2]*c(-1,1),
				  lwd=lwdV,col=colV)
		}
	}
}

VisualiseXYV_local <- function(X,Y=NULL,V=NULL,N,magfac=0.5,
		x0=c(median(X[,1]),median(X[,2])),
		ltyx0=1,colx0='red',
		xsf=1,offset=0.07)
	# Input are data matrices
	# -  X (n x 2) : containing locations X[i,],
	# -  Y (n x 2) : containing direction vectors Y[i,] at X[i,]
	#                (not necessarily unit vectors),
	# -  V (n x 2) : containing axis vectors V[i,] (and
	#                implicitly -V[i,]) at X[i,]
	#                (not necessarily unit vectors).
	# The procedure plots the locations XX[i,] and the
	# corresponding directions Y[i,] and axes V[i,]
	# on a gray scale reflecting the weights W[i] of the
	# observations as computed by SetWeights(S,N), where
	# S[i] == sum((X[i,] - x0)^2).
	# The reference point x0 can be changed interactively
	# by clicking into the plot window. (To terminate, one
	# has to click outside.)
{
	n <- dim(X)[1]
	if (N >= n | N <= 0){
		return('Input parameter N makes no sense!')
	}
	par(cex=1.2,mai=c(0.7,0.75,0.1,0.1),mgp=c(1.7,0.5,0))
	# 1st plot:
	S <- colSums(((t(X) - x0)/xsf)^2)
	W <- SetWeights(S,N)
	colX=gray((1-offset/2)*(1-W))
	colYV=gray((1-offset)*(1-W))
	plot(X[,1],X[,2],pch=16,col=colX,
		 xlab=expression(italic(x[1])),
		 ylab=expression(italic(x[2])),
		 main='')
	if (!is.null(Y)){
		for (i in (1:n)){
			lines(X[i,1] + magfac*Y[i,1]*c(0,1),
				  X[i,2] + magfac*Y[i,2]*c(0,1),
				  col=colYV[i])
		}
	}
	if (!is.null(V)){
		for (i in (1:n)){
			lines(X[i,1] + magfac*V[i,1]*c(-1,1),
				  X[i,2] + magfac*V[i,2]*c(-1,1),
				  col=colYV[i])
		}
	}
	abline(v=x0[1],h=x0[2],lty=ltyx0,col=colx0)
	# New point:
	x0old <- x0
	x0 <- locator(1)
	x0 <- c(x0$x,x0$y)
	while (x0[1] >= min(X[,1]) & x0[1] <= max(X[,1]) &
		   x0[2] >= min(X[,2]) & x0[2] <= max(X[,2])){
		# 1st plot:
		S <- colSums(((t(X) - x0)/xsf)^2)
		W <- SetWeights(S,N)
		colX=gray((1-offset/2)*(1-W))
		colYV=gray((1-offset)*(1-W))
		plot(X[,1],X[,2],pch=16,col=colX,
			 xlab=expression(italic(x[1])),
			 ylab=expression(italic(x[2])),
			 main='')
		if (!is.null(Y)){
			for (i in (1:n)){
				lines(X[i,1] + magfac*Y[i,1]*c(0,1),
					  X[i,2] + magfac*Y[i,2]*c(0,1),
					  col=colYV[i])
			}
		}
		if (!is.null(V)){
			for (i in (1:n)){
				lines(X[i,1] + magfac*V[i,1]*c(-1,1),
					  X[i,2] + magfac*V[i,2]*c(-1,1),
					  col=colYV[i])
			}
		}
		abline(v=x0[1],h=x0[2],lty=ltyx0,col=colx0)
		# New point:
		x0old <- x0
		x0 <- locator(1)
		x0 <- c(x0$x,x0$y)
	}
	return(x0old)
}


# Auxiliary function ----

GreedySubset <- function(X,size,V=NULL, returnJ=FALSE)
	# Auxiliary function.
	# For a data matrix X (n x q) with n,q > 1, this procedure
	# determines stepwise a submatrix X0 = X[J,] with size rows
	# such that the single rows are approximately as far apart
	# as possible.
{
	n <- dim(X)[1]
	if (n <= size){
		if (is.null(V)){
			return(X)
		}else{
			return(list(X0=X,V0=V))
		}
	}
	J <- rep(NA,size)
	J[1] <- sample(n,1)
	dist <- colSums((t(X) - X[J[1],])^2)
	for (i in 2:size){
		J[i] <- which.max(dist)
		dist <- pmin(dist,colSums((t(X) - X[J[i],])^2))
	}
	J <- sort(J)
	if (returnJ==FALSE){
		if (is.null(V)){
			return(X[J,])
		}else{
		return(list(X0=X[J,],V0=V[J,]))
		}
	}else{
		# we also return J
		if (is.null(V)){
			return(list(X[J,], J=J))
		}else{
		return(list(X0=X[J,],V0=V[J,], J=J))
		}
	}
}
