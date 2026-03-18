# Auxiliary function G and its first two derivatives ----

# G(t) = sum_{k >= 0} c(d,k) t^{2k}
# with c(d,0) = 1 and
# c(d,k+1)/c(d,k) = 1/[(2k + 2)(2k + d).

# In case of d = 3,
# G(t) = sinh(t)/t.

G0 <- function(t,d=3,prec=10^(-11))
{
	k <- 0
	At <- 1
	G0t <- 1
	t24 <- t^2/4
	Ft <- t24/(k+1)/(k+d/2)
	while (max(Ft) > 0.5 || max(At) > prec){
		k <- k+1
		At <- At*Ft
		G0t <- G0t + At
		Ft <- t24/(k+1)/(k+d/2)
	}
	return(G0t)
}

G1 <- function(t,d=3,prec=10^(-11))
{
	G1t <- (t/d)*G0(t,d+2,prec)
	return(G1t)
}

G2 <- function(t,d=3,prec=10^(-11))
{
	G2t <- G0(t,d+2,prec)/d + (t^2/d/(d+2))*G0(t,d+4,prec)
	return(G2t)
}

G01 <- function(t,d=3,prec=10^(-11))
{
	G0t <- G0(t,d,prec)
	G1t <- (t/d)*G0(t,d+2,prec)
	return(list(G0t=G0t,G1t=G1t))
}

G012 <- function(t,d=3,prec=10^(-11))
{
	G0t <- G0(t,d,prec)
	G0tmp <- G0(t,d+2,prec)
	G1t <- (t/d)*G0tmp
	G2t <- G0tmp/d + (t^2/d/(d+2))*G0(t,d+4,prec)
	return(list(G0t=G0t,G1t=G1t,G2t=G2t))
}


lG0 <- function(t,d=3,prec=10^(-11),threshold=100*d)
{
	lgt <- rep(0,length(t))
	tmp <- (t < threshold)
	if (any(tmp)){
		lgt[tmp] <- log(G0(t[tmp],d,prec))
	}
	tmp <- !tmp
	if (any(tmp)){
		ad <- (d-1)/2
		Oo2ttmp <- 1/2/t[tmp]
		lgt[tmp] <- (ad - 1)*log(2) +
			lgamma(d/2) - lgamma(1/2) +
			t[tmp] - ad*log(t[tmp]) -
			ad*(ad-1)*Oo2ttmp*(1 + Oo2ttmp)
	}
	return(lgt)
}

# Moments of the vMF distribution ----

# We parametrize the Von Mises-Fisher distribution on the unit sphere of R^d by a vector z in R^d.

Moments1_vMF <- function(t,d=3,prec=10^(-11),threshold=100*d)
	# For t >= 0, this procedure computes the quantities
	#    M0 = log(G_d(t)),
	#    M1 = E(U_t),
	#    M2a = Var(U_t),
	#    M2b = (1 - E(U_t)^2)/(d-1)
	# where U_t has density
	#    exp(tu - tilde{gamma}_d(t))*h_d(u).
	# The latter three quantitites are needed for the expectation
	# and covariance of Q_z, where z is a vector in R^d with
	# norm t:
	#    mean(Q_z) = M1 * v ,
	#    cov(Q_z)  = M2a * P + M2b * (I_d - P)
	# with P = v %*% t(v) and v = z/t.
{
	M0 <- rep(0,length(t))
	M1 <- rep(0,length(t))
	M2a <- rep(0,length(t))
	M2b <- rep(0,length(t))
	# Exact formulae for moderate values of t:
	tmp <- (t < threshold)
	if (any(tmp)){
		G012ttmp <- G012(t[tmp],d,prec)
		M0[tmp] <- log(G012ttmp$G0t)
		M1[tmp] <- G012ttmp$G1t/G012ttmp$G0t
		M2ttmp <- G012ttmp$G2t/G012ttmp$G0t
		M2a[tmp] <- M2ttmp - M1[tmp]^2
		M2b[tmp] <- (1 - M2ttmp)/(d-1)
	}
	# Approximate formulae for large values of t:
	tmp <- !tmp
	if (any(tmp)){
		ad <- (d-1)/2
		Oo2ttmp <- 1/2/t[tmp]
		M0[tmp] <- (ad - 1)*log(2) +
			lgamma(d/2) - lgamma(1/2) +
			t[tmp] - ad*log(t[tmp]) -
			ad*(ad-1)*Oo2ttmp*(1 + Oo2ttmp)
		M1[tmp] <- 1 - ad/t[tmp]*(1 - (ad-1)*Oo2ttmp)
		M2a[tmp] <- ad/t[tmp]^2*(1 - (ad-1)/t[tmp])
		M2b[tmp] <- 1/t[tmp]*(1 - ad/t[tmp])
	}
	
	return(list(M0=M0,M1=M1,M2a=M2a,M2b=M2b))
}

Moments_vMF <- function(z,d=length(z),prec=10^(-11),
						threshold=100*d)
	# This procedure computes a list with three items:
	#  gammaz = gamma_d(z) = log(G_d(||z||)),
	#  muz = mean(Q_z),
	#  Sigmaz = cov(Q_z).
{
	if (length(z) != d){
		print('Oops, dimension of z is different from d!?')
		return()
	}
	t <- sqrt(sum(z^2))
	if (t==0){
		return(list(gammaz=0,muz=rep(0,d),Sigmaz=diag(d)/d))
	}
	v <- z/t
	Mt <- Moments1_vMF(t,d,prec,threshold)
	gammaz <- Mt$M0
	muz <- Mt$M1 * v
	Sigmaz <- Mt$M2a * (v %*% t(v)) +
			Mt$M2b * (diag(d) - v %*% t(v))
	return(list(gammaz=gammaz,muz=muz,Sigmaz=Sigmaz))
}

Mu_vMF <- function(z,d=length(z),prec=10^(-11),threshold=100*d)
{
	if (length(z) != d){
		print('Oops, dimension of z is different from d!?')
		return()
	}
	t <- sqrt(sum(z^2))
	if (t==0){
		return(z)
	}
	v <- z/t
	if (t < threshold){
		G01t <- G01(t,d,prec)
		muz <- (G01t$G1t / G01t$G0t) * v
	}else{
		ad <- (d-1)/2
		muz <- (1 - ad/t*(1 - (ad-1)/2/t))*v
	}
	return(muz)
}


# Maximum-likelihood estimation ----

# If (w_1,Y_1), ..., (w_n,Y_n) are independent observations with Y_i ~ vMF(z) for some unknown parameter z in R^d and fixed weights w_1, ..., w_n > 0, then the unique minimizer of the weighted negative log-likelihood
#   L(z) = sum of w_i (gamma(z) - Y_i'z)
# is the unique vector z such that
#   mu(z) = Yw = sum_i w_i Y_i / sum_i w_i.
# Here gamma(z) = log(H0(||z||^2)).
# It can be computed by calling MLE_vMF(Yw,d).

MLE0_vMF <- function(y,d=length(y),prec=10^(-11),
					 rmax=1-(1 - 1/d)/200)
# Writing y = ru with a unit vector u and a scalar r >= 0,
# this procedure compute the solution z = tu of
#    Mu_vMF(z,d,prec) = y
# by minimizing log(G_d(t)) - rt.
# For r >= rmax, an approximation is used.
{
	r <- sqrt(sum(y^2))
	if (r >= 1){
		print('There exists no solution!')
		return()
	}
	if (r == 0){
		return(y)
	}
	if (r > rmax){
		ad <- (d-1)/2
		t <- ad/2/(1-r)
		t <- t*(1 + sqrt(1 - (ad-1)/t))
		return(t*(y/r))
	}else{
		t <- d*r
		Mt <- Moments1_vMF(t,d,prec,Inf)
		F1t <- Mt$M1 - r
		iter1 <- 0
		while (abs(F1t) > prec & iter1 < 1000){
			iter1 <- iter1 + 1
			F2t <- Mt$M2a
			tnew <- t - F1t/F2t
			F0t <- Mt$M0 - r*t
			F0tnew <- lG0(tnew,d,prec,Inf) - r*tnew
			iter0 <- 0
			while (iter0 < 20 & F0tnew >= F0t){
				iter0 <- iter0 + 1
				tnew <- (t + tnew)/2
				s <- tnew^2
				F0tnew <- lG0(tnew,d,prec,Inf) - r*tnew
			}
			if (iter0 == 20){
				F1t <- 0
			}else{
				t <- tnew
				Mt <- Moments1_vMF(t,d,prec,Inf)
				F1t <- Mt$M1 - r
			}
		}
		return(t*(y/r))
	}
}


MLE_vMF <- function(y,d=length(y),prec=10^(-11),threshold=100*d)
	# Minimize L(z) = log(G_d(||z||) - y'z.
	# This procedure is not very efficient, but it may be
	# generalized for regression settings.
{
	if (sum(y^2) >= 1){
		print('There exists no solution!')
		return()
	}
	z <- d*y
	Mz <- Moments_vMF(z,d,prec,threshold)
	Lz <- Mz$gammaz - sum(y*z)
	muz <- Mz$muz
	Sigmaz <- Mz$Sigmaz
	dz <- qr.solve(Sigmaz,y-muz)
	dirderiv <- sum(dz*(y-muz))
	iter1 <- 0
	while (dirderiv > prec & iter1 < 1000){
		iter1 <- iter1 + 1
		znew <- z + dz
		tnew <- sqrt(sum(znew^2))
		Lznew <- lG0(tnew,d,prec,threshold) - sum(y*znew)
		iter0 <- 0
		while (iter0 < 20 & Lznew >= Lz){
			iter0 <- iter0 + 1
			znew <- (z + znew)/2
			tnew <- sqrt(sum(znew^2))
			Lznew <- lG0(tnew,d,prec,threshold) - sum(y*znew)
		}
		if (iter0 == 20){
			dirderiv <- 0
		}else{
			z <- znew
			Mz <- Moments_vMF(z,d,prec,threshold)
			Lz <- Mz$gammaz - sum(y*z)
			muz <- Mz$muz
			Sigmaz <- Mz$Sigmaz
			dz <- qr.solve(Sigmaz,y-muz)
			dirderiv <- sum(dz*(y-muz))
		}
	}
	return(z)
}

# Simulate samples from vMF distributions ----

Sample_vMF <- function(n,z,d=length(z))
{
	mz <- sqrt(sum(z^2))
	Y <- matrix(0,n,d)
	n1 <- 0
	n2 <- n
	while (n2 > 0){
		if (n2 == 1){
			Ynew <- rnorm(d)
			Ynew <- Ynew/sqrt(sum(Ynew^2))
			if (runif(1) <= exp(sum(Ynew*z)-mz)){
				Y[n,] <- Ynew
				n2 <- 0
			}
		}else{
			Ynew <- matrix(rnorm(n2*d),n2,d)
			Ynew <- Ynew/sqrt(rowSums(Ynew^2))
			done <- (runif(n2) <= exp(colSums(t(Ynew)*z) - mz))
			n0 <- sum(done)
			if (n0 > 0){
				Y[(n1+1):(n1+n0),] <- Ynew[done,]
				n1 <- n1+n0
				n2 <- n2-n0
			}
		}
	}
	return(Y)
}
