# Nonparametric isotonic fitting ----

IsotonicH <- function(AJ)
	# Input is a data matrix AJ with N rows and 2m columns
	# A1, J1, A2, J2, ..., Am, Jm,
	# where Ak[i] is an angle in [0, pi/2]
	# and Jk[i] is a number in [0,1].
	# The program computes a non-decreasing function
	# h : [0,pi/2] --> [0,1]
	# such that
	#    sum_{i=1}^N sum_{k=1}^m (Jk[i] - h(Ak[i]))^2
	# is minimal. The output is a matrix TH with
	# TH[,1] containing the elements t of
	# {Ak[i] : i <= N, k <= m} in increasing order
	# and TH[,2] containing the corresponding values h(t).
{
	N <- dim(AJ)[1]
	m <- dim(AJ)[2]/2
	A <- AJ[,seq(1,2*m-1,2)]
	J <- AJ[,seq(2,2*m,2)]
	t <- sort(unique(A))
	n <- length(t)
	w <- rep(0,n)
	y <- rep(0,n)
	for (i in 1:n){
		tmp <- (A == t[i])
		w[i] <- sum(tmp)
		y[i] <- mean(J[tmp])
	}
	ht <- IsoMeans(y,w)
	return(cbind('t'=t,'h(t)'=ht))
}

# Auxiliary function (PAVA)

IsoMeans <- function(y, w=rep(1,times=length(y)))
	# IsoMeans(y) fits an isotonic vector mhat to a data
	# vector y such that
	#    sum((y - mhat)^2) is minimal
	# (Pool-adjacent-violators algorithm).
	# In case of a second input vector w with positive
	# entries (and the same size as y), IsoMeans(y,w)
	# produces an isotonic vector mhat minimizing
	# sum(w*(y - mhat)^2).
	# 
	# Lutz Duembgen, October 3, 2008
{
	
	n <- length(y)
	index <- rep(0,times=n)
	weight <- rep(0,times=n)
	# An interval of indices is represented by its left
	# endpoint ("index") and its "weight" (sum of w's).
	mhat <- rep(0,times=n)
	
	ci <- 1
	index[ci] <- 1
	j <- 1
	# ci is the number of the interval considered currently,
	# where we start with {1}. The index j is the right endpoint
	# of the current interval.
	weight[ci] <- w[1]
	mhat[ci] <- y[1]
	# mhat[ci] is the weighted mean of y-values within
	# the current interval.
	while (j < n)
	{
		j <- j + 1
		# a new index intervall {j} is created:
		ci <- ci + 1
		index[ci] <- j
		weight[ci] <- w[j]
		mhat[ci] <- y[j]
		while (ci >= 2 && mhat[ci] <= mhat[ci-1])
		{
			# "pool adjacent violators":
			nw <- weight[ci-1] + weight[ci]
			mhat[ci-1] <- mhat[ci-1] +
				(weight[ci] / nw) * (mhat[ci] - mhat[ci-1])
			weight[ci-1] <- nw
			ci <- ci-1
		}
	}
	# Now define mhat for all indices:
	while (j >= 1)
	{
		mhat[index[ci]:j] <- mhat[ci]
		j <- index[ci] - 1
		ci <- ci - 1
	}
	return(mhat)
}

# Parametric fitting with (natural) cubic splines ----

FittedH <- function(AJ,xx = seq(0,pi/2,length.out=181),
					M=4,end.at.1=TRUE,trunc.at.1=TRUE)
	# Input is a data matrix AJ with N rows and 2m columns
	# A1, J1, A2, J2, ..., Am, Jm,
	# where Ak[i] is an angle in [0, pi/2]
	# and Jk[i] is a number in [0,1].
	# The program computes a cubic spline function h on
	# [0,pi/2] with M+1 equidistant knots partitioning in
	# this interval such that h' = 0 at the boundaries and
	#    sum_{i=1}^N sum_{k=1}^m (Jk[i] - h(Ak[i]))^2
	# is minimimal.
	# If end.at.1 == TRUE, then the additional constraint
	# h(pi/2) = 1 is used.
	# If trunc.at.1 == TRUE, then function values > 1 are
	# replaced with 1.
	# The output is a matrix XH with
	# XH[,1] containing the given vector xx and
	# XH[,2] the corresponding values h(xx[i]).
{
	N <- dim(AJ)[1]
	m <- dim(AJ)[2]/2
	A <- as.vector(AJ[,seq(1,2*m-1,2)])
	J <- as.vector(AJ[,seq(2,2*m,2)])
	BBA <- CubicSplines0(M=M,xx=A)
	BBxx <- CubicSplines0(M=M,xx=xx)
	if (end.at.1){
		BBA1 <- BBA
		BBA1[,M] <- BBA1[,M] - 0.5*BBA1[,M+1]
		BBA1[,M+1] <- 1.5*BBA1[,M+1]
		J1 <- J - BBA1[,M+1]
		theta <- c(qr.solve(BBA1[,1:M],J1),1)
		BBxx1 <- BBxx
		BBxx1[,M] <- BBxx1[,M] - 0.5*BBxx1[,M+1]
		BBxx1[,M+1] <- 1.5*BBxx1[,M+1]
		hxx <- BBxx1 %*% theta
	}else{
	theta <- qr.solve(BBA,J)
	hxx <- BBxx %*% theta
	}
	if (trunc.at.1){
		hxx <- pmin(hxx,1)
	}
	return(cbind('x'=xx,'h(x)'=hxx))
}

# Auxiliary functions:

BSplines <- function(t,d=3,
					 xx=seq(t[1],t[length(t)],length.out=401),
					 delta.factor=1)
	# Returns a matrix BB with n=length(xx) rows and
	# m+d columns (m = length(t) -1), where the j-th
	# column of BB contains the j-th basis function
	# evaluated at xx.
{
	m <- length(t) - 1
	n <- length(xx)
	
	# Augmented vector of knots:
	tmp <- median(t[2:(m+1)] - t[1:m])*delta.factor
	t <- c(t[1] - (d:1)*tmp, t, t[m+1] + (1:d)*tmp)
	
	# Start with the B-spline basis of order 0:
	mm <- m + 2*d
	BB <- matrix(0,n,mm)
	for (j in 1:mm)
	{
		BB[xx >= t[j] & xx < t[j+1], j] <- 1
	}
	
	# Compute inductively higher order bases:
	for (j in 1:d)
	{
		CC <- matrix(0,n,mm-j)
		for (z in 1:(mm-j))
		{
			d0 <- t[z+j] - t[z]
			d1 <- t[z+1+j] - t[z+1]
			if (d0 > 0)
			{
				CC[,z] <- CC[,z] +
					(xx - t[z])/d0 * BB[,z]
			}
			if (d1 > 0)
			{
				CC[,z] <- CC[,z] +
					(t[z+1+j] - xx)/d1 * BB[,z+1]
			}
		}
		BB <- CC
	}
	
	if (delta.factor == 0)
	{
		BB[n,m+d] <- 1
	}
	
	return(BB)
}

CubicSplines0 <- function(tmin=0,tmax=pi/2,M=4,
					 xx=seq(tmin,tmax,length.out=181))
	# Cubic splines on [tmin,tmax] with M+1 equidistant
	# knots and derivative 0 at the boundaries.
	# Returns a matrix BB with n=length(xx) rows and
	# M+1 columns, where the j-th column of BB contains
	# the j-th basis function evaluated at xx.
{
	t <- seq(tmin,tmax,length.out=M+1)
	BB <- BSplines(t,d=3,xx=xx,delta.factor=1)
	BB[,3] <- BB[,1] + BB[,3]
	m <- length(t) - 1
	BB[,m+1] <- BB[,m+1] + BB[,m+3]
	BB <- BB[,2:(m+2)]
	return(BB)
}


# CH: hardcoded function to read in the file and produce the spline fit
GetSunFit <- function(dcpath, M=3, Jshifta=NULL, Jshiftb=NULL, Jshiftc=NULL){
	# M is passed to FittedH
	# the scaleJ parameter can be used to scale J, currently with numbers 0, 0.33, 0.67, 1
	# note, file is here: setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
	dc <- read.csv(dcpath)
	# bring dc rows into matrix form
	A <- c(dc$A00, dc$A45, dc$A90, dc$A135)
	J <- c(dc$J00, dc$J45, dc$J90, dc$J135)
	# shift if desired:
	if(!is.null(Jshifta)){
		# former 0
		J[J==0] <- Jshifta
	}
	if(!is.null(Jshiftb)){
		# former 0.33
		J[J==0.33] <- Jshiftb
	}
	if(!is.null(Jshiftc)){
		# former 0.67
		J[J==0.67] <- Jshiftc
	}
	# NOTE: no shfit for 1

	# modify so that Aik element [0, pi/2] instead of [0, pi]
	A[A>pi/2] <- pi - A[A>pi/2]
	N <- length(dc$J45)
	m <- 4
	AJ <- matrix(0,N,2*m) # 2*m, 1 for A, 1 for J
	AJ[,seq(1,2*m-1,2)] <- A
	AJ[,seq(2,2*m,2)] <- J
	# fit the data:
	XH <- FittedH(AJ,M=M)
	# the weights are 1/h(t)
	return(cbind(XH[,1],1/XH[,2]))
}