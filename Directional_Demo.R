source("Directional.R")

# Moments of the vMF distribution ----

# We parametrize the Von Mises-Fisher distribution on the unit sphere of R^d by a vector z in R^d.

z <- rnorm(3,0,0.1)
z
d <- length(z)
Mu_vMF(z)
cbind(z,Mu_vMF(z),1/d*(1 - sum(z^2)/d/(d+2))*z)

z <- rnorm(3,0,2)
d <- length(z)
t <- sqrt(sum(z^2))
M0 <- Mu_vMF(z)
s0 <- 0.00001
lG1 <- lG0(t,d)
lG2 <- rep(0,d)
for (j in 1:length(z)){
	zj <- z
	zj[j] <- zj[j] + s0
	lG2[j] <- lG0(sqrt(sum(zj^2)),d)
}
cbind(M0, (lG2 - lG1)/s0)

S1 <- Moments_vMF(z)$Sigmaz
eigen(S1)$values

s0 <- 0.000001
M0 <- Mu_vMF(z)
M <- matrix(0,d,d)
for (j in 1:length(z)){
	zj <- z
	zj[j] <- zj[j] + s0
	M[,j] <- Mu_vMF(zj)
}
S2 <- (M - M0)/s0
S2 - S1

# Maximum-likelihood estimation ----

# If (w_1,Y_1), ..., (w_n,Y_n) are independent observations with Y_i ~ vMF(z) for some unknown parameter z in R^d and fixed weights w_1, ..., w_n > 0, then the unique minimizer of the weighted negative log-likelihood
#   L(z) = sum of w_i (gamma(z) - Y_i'z)
# is the unique vector z such that
#   mu(z) = Ybar = sum_i w_i Y_i / sum_i w_i.
# Here gamma(z) = log(H0(||z||^2,d)).
# It can be computed by calling MLE_vMF(Yw).

z <- c(7,1,0)
n <- 200
Y <- Sample_vMF(n,z)
plot(cos(2*pi*seq(0,1,0.001)),sin(2*pi*seq(0,1,0.001)),
	 type='l',col='gray',
	 xlab='Y[1]',ylab='Y[2]')
# If d = 2:
points(Y[,1],Y[,2],pch=16)
# If d >= 3:
oY <- order(Y[,3])
points(Y[oY,1],Y[oY,2],pch=19,col=gray((1 - Y[oY,3])/2))
Ybar <- colMeans(Y)
Ybar
sqrt(sum(Ybar^2))

d <- 3
MLE0_vMF(c(0.998,rep(0,d-1)))
MLE_vMF(c(0.998,rep(0,d-1)))

MLE0_vMF(Ybar)
MLE_vMF(Ybar)
Mu_vMF(MLE0_vMF(Ybar)) - Ybar
Mu_vMF(MLE_vMF(Ybar)) - Ybar
