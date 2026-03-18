setwd('C:/Users/chaslebacher/Desktop/R_code/Europa_Bingham')

source('Directional.R')
source('GLM_3D.R')
source('AxialDataSphere.R')

# Illustrate stereographic projections ----

set.seed(3032024)
n1 <- 20
n2 <- 100
n <- n1*n2
X <- matrix(0,n,3)
V <- matrix(0,n,3)
for (i in 1:n1){
	II <- ((i-1)*n2+1):(i*n2)
	u <- runif(1,-1,1)
	u2 <- sqrt(1 - u^2)
	dir <- rnorm(3)
	dir <- dir/sqrt(sum(dir^2))
	Xi <- matrix(rnorm(n2*3),n2,3)
	Xi <- t(t(Xi) - dir %*% t(Xi %*% dir))
	Xi <- Xi/sqrt(rowSums(Xi^2))
	Xi <- t(u2*t(Xi) + u*dir)
	X[II,] <- Xi
	V[II,] <- t(dir - u*t(Xi)) / sqrt(1 - u^2)
}
str(X)
str(V)
range(rowSums(X^2))
range(rowSums(V^2))
range(rowSums(X*V))

par(cex=1.2,mai=c(0.5,0.5,0.1,0.1))
PlotAxialDataSphere(X=X,V=V)
x0 <- ShowAxialDataSphere(X,V,vfactor=0.05)
res <- StereogrProj(x0=x0,X=X,V=V,xysize=1,vfactor=0.05)
res <- StereogrProj(x0=x0,X=X,V=V,xysize=2,vfactor=0.025)
res <- StereogrProj(x0=x0,X=X,V=V,xysize=4,vfactor=0.0125)
par(cex=1.2,mai=c(0.5,0.5,0.5,0.1))
x0 <- StereogrTour(X=X,V=V,x0=x0)

str(res)
range(rowSums(res$V2^2)-1)

X2 <- res$X2
V2 <- res$V2

res2 <- InvStereogrProj(x0=x0,X2=X2,V2=V2)
str(res2)

range(rowSums(X * V))
range(rowSums(res2$X * res2$V))
range(rowSums(res2$X^2) - 1)
range(rowSums(res2$V^2) - 1)

range(X - res2$X)
range(V - res2$V)

# Analyze data from Europa ----
ds <- read.csv("ULDR_children_for_vMF_spherical.csv",
			   header=TRUE,sep=',')
str(ds)

X <- as.matrix(ds[,2:4])
str(X)
V <- as.matrix(ds[,8:10])
str(V)
range(rowSums(X^2))
range(rowSums(V^2))
range(rowSums(X*V))

C <- ds$cluster
str(C)
table(C)

subsetsize <- 200
n0 <- sum(pmin(table(C),subsetsize))
n0

# Illustrate local weights only for cluster C==7:
x0 <- ShowWeightedAxialDataSphere(X[C==7,],V[C==7,],N=200,
								  vfactor=0.05)
x0
x0 <- ShowWeightedAxialDataSphere(X[C==7,],V[C==7,],N=200,
								  vfactor=0.05,x0=x0)
points(0,0,pch=19,col='red')

X0 <- matrix(0,n0,3)
V0 <- matrix(0,n0,3)
W0const <- matrix(0,n0,3)
M0const <- matrix(0,n0,3)
W0lin   <- matrix(0,n0,3)
M0lin   <- matrix(0,n0,3)
W0quadr <- matrix(0,n0,3)
M0quadr <- matrix(0,n0,3)

# N <- 400
# a0 <- 0
# set.seed(3032024)
# for (cl in 0:18){
# 	tmp <- which(C==cl)
# 	tmpXV <- GreedySubset(X=X[tmp,],size=subsetsize,V=V[tmp,])
# 	print(paste('Cluster',cl))
# 	b0 <- dim(tmpXV$X0)[1]
# 	X0[a0 + (1:b0),] <- tmpXV$X0
# 	V0[a0 + (1:b0),] <- tmpXV$V0
# 	res <- SmoothAxialDataSphere(X0=tmpXV$X0,
# 	 							 X=X[tmp,],V=V[tmp,],
#  	 							 N=N)
# 	W0const[a0 + (1:b0),] <- res$W0const
# 	M0const[a0 + (1:b0),] <- res$M0const
# 	W0lin[a0 + (1:b0),]   <- res$W0lin
# 	M0lin[a0 + (1:b0),]   <- res$M0lin
# 	W0quadr[a0 + (1:b0),] <- res$W0quadr
# 	M0quadr[a0 + (1:b0),] <- res$M0quadr
# 	a0 <- a0 + b0
# }

# # Save results:
# Europa <- data.frame('X0'=X0,'V0'=V0,
# 					'W0const'=W0const,'M0const'=M0const,
# 						'W0lin'=W0lin,'M0lin'=M0lin,
# 						'W0quadr'=W0quadr,'M0quadr'=M0quadr)
# saveRDS(Europa,file='Europa400.rds')

# Recover results:
Europa <- readRDS(file='Europa100.rds')
Europa <- readRDS(file='Europa200.rds')
Europa <- readRDS(file='Europa400.rds')
str(Europa)
X0 <- as.matrix(Europa[,1:3])
V0 <- as.matrix(Europa[,4:6])
W0const <- as.matrix(Europa[,7:9])
M0const <- as.matrix(Europa[,10:12])
W0lin   <- as.matrix(Europa[,13:15])
M0lin   <- as.matrix(Europa[,16:18])
W0quadr <- as.matrix(Europa[,19:21])
M0quadr <- as.matrix(Europa[,22:24])

x0 <- ShowWeightedAxialDataSphere(X,V,N=100)
x0 <- ShowWeightedAxialDataSphere(X,V,N=200,
								  vfactor=0.05,x0=x0)


par(cex=1.2,mai=c(0.5,0.5,0.1,0.1))
# Look at the data and pick an interesting point x0:
# Complete raw data
x0 <- ShowAxialDataSphere(X=X,V=V,vfactor=0.02)
# Subset of raw data:
x0 <- ShowAxialDataSphere(X=X0,V=V0,vfactor=0.05)
# Only locations:
x0 <- ShowAxialDataSphere(X=X,V=NULL)
# Only subset of locations:
x0 <- ShowAxialDataSphere(X=X0,V=NULL)
# For paper: 10 x right, 3 x up.
# Local constant fits:
x0 <- ShowAxialDataSphere(X=X0,V=M0const,x0=x0,
						  Vcolours = c('blue','gray'))
# Local linear fits:
x0 <- ShowAxialDataSphere(X=X0,V=M0lin,x0=x0,
						  Vcolours = c('blue','gray'))
# Local quadratic fits:
x0 <- ShowAxialDataSphere(X=X0,V=M0quadr,x0=x0,
						  Vcolours = c('blue','gray'))

# Pick a size for the projection window:
xysize <- 0.45
# Look at the stereographic projection:
# Complete data:
tmp <- StereogrProj(x0=x0,X=X,V=V,xysize=xysize)
# Only subset of data:
tmp <- StereogrProj(x0=x0,X=X0,V=V0,xysize=xysize) # subset
# Only subset of locations without axes:
tmp <- StereogrProj(x0=x0,X=X0,xysize=xysize)
# Local constant fits:
tmp <- StereogrProj(x0=x0,X=X0,V=M0const,xysize=xysize)
# Local linear fits:
tmp <- StereogrProj(x0=x0,X=X0,V=M0lin,xysize=xysize)
# Local quadratic fits:
tmp <- StereogrProj(x0=x0,X=X0,V=M0quadr,xysize=xysize)


# Analyze data from Ganymede, V1 ----

ds <- read.csv("Ganymede_CRossi_for_vMF_spherical.csv",
			   header=TRUE,sep=',')
str(ds)

X <- as.matrix(ds[,2:4])
str(X)
V <- as.matrix(ds[,8:10])
str(V)
range(rowSums(X^2))
range(rowSums(V^2))
range(rowSums(X*V))

# Look at the data:
x0 <- ShowAxialDataSphere(X=X,V=V,dangle=2*pi/24)

# tmp <- GreedySubset(X,size=2000,V=V)
# str(tmp)
# X0 <- tmp$X0
# V0 <- tmp$V0
# str(X0)
# str(V0)
# range(rowSums(X0^2))
# range(rowSums(V0^2))
# range(rowSums(X0*V0))
# 
# # Look at a subset of the data:
# x0 <- ShowAxialDataSphere(X=X0,V=V0)
# 
# PF <- SmoothAxialDataSphere(X0=X0,X=X,V=V,N=200)
# str(PF)
# 
# # Save results:
# Ganymede400 <- data.frame('X0'=X0,'V0'=V0,
# 	'W0const'=PF$W0const,'M0const'=PF$M0const,
# 	'W0lin'=PF$W0lin,'M0lin'=PF$M0lin,
# 	'W0quadr'=PF$W0quadr,'M0quadr'=PF$M0quadr)
# saveRDS(Ganymede400,file='Ganymede400.rds')

# Recover results:
Ganymede <- readRDS(file='Ganymede100.rds')
Ganymede <- readRDS(file='Ganymede200.rds')
Ganymede <- readRDS(file='Ganymede400.rds')
str(Ganymede)

X0 <- as.matrix(Ganymede[,1:3])
V0 <- as.matrix(Ganymede[,4:6])
W0const <- as.matrix(Ganymede[,7:9])
M0const <- as.matrix(Ganymede[,10:12])
W0lin   <- as.matrix(Ganymede[,13:15])
M0lin   <- as.matrix(Ganymede[,16:18])
W0quadr <- as.matrix(Ganymede[,19:21])
M0quadr <- as.matrix(Ganymede[,22:24])

# Show subset of data:
x0 <- ShowAxialDataSphere(X=X0,V=V0)
# Show local constant fits:
x0 <- ShowAxialDataSphere(X=X0,V=M0const)
# Show local linear fits:
x0 <- ShowAxialDataSphere(X=X0,V=M0lin)
# Show local quadratic fits:
x0 <- ShowAxialDataSphere(X=X0,V=M0quadr)

# Analyze data from Ganymede, V2 ----

ds <- read.csv("Ganymede_CRossi_for_vMF_spherical.csv",
			   header=TRUE,sep=',')
str(ds)

X <- as.matrix(ds[,2:4])
str(X)
V <- as.matrix(ds[,8:10])
str(V)
range(rowSums(X^2))
range(rowSums(V^2))
range(rowSums(X*V))

X0 <- read.csv("Ganymede_grid_for_R.csv",
			   header=TRUE,sep=',')
str(X0)
X0 <- as.matrix(X0)
str(X0)

# PF <- SmoothAxialDataSphere(X0=X0,X=X,V=V,N=100)
# PF <- SmoothAxialDataSphere(X0=X0,X=X,V=V,N=200)
# PF <- SmoothAxialDataSphere(X0=X0,X=X,V=V,N=400)
# str(PF)
# 
# range(rowSums(X0*PF$W0const))
# range(rowSums(X0*PF$M0const))
# range(rowSums(X0*PF$W0lin))
# range(rowSums(X0*PF$M0lin))
# range(rowSums(X0*PF$W0quadr))
# range(rowSums(X0*PF$M0quadr))
# 
# x0 <- ShowAxialDataSphere(X=X0,V=PF$M0const,vfactor=0.1)
# x0 <- ShowAxialDataSphere(X=X0,V=PF$M0lin,vfactor=0.1)
# x0 <- ShowAxialDataSphere(X=X0,V=PF$M0quadr,vfactor=0.1)

# # Save results:
# Ganymede <- data.frame('X0'=X0,
# 	'W0const'=PF$W0const,'M0const'=PF$M0const,
# 	'W0lin'=PF$W0lin,'M0lin'=PF$M0lin,
# 	'W0quadr'=PF$W0quadr,'M0quadr'=PF$M0quadr)
# saveRDS(Ganymede,file='GanymedeB100.rds')
# saveRDS(Ganymede,file='GanymedeB200.rds')
# saveRDS(Ganymede,file='GanymedeB400.rds')

# Recover results:
Ganymede <- readRDS(file='GanymedeB100.rds')
Ganymede <- readRDS(file='GanymedeB200.rds')
Ganymede <- readRDS(file='GanymedeB400.rds')
str(Ganymede)

X0 <- as.matrix(Ganymede[,1:3])
W0const <- as.matrix(Ganymede[,4:6])
M0const <- as.matrix(Ganymede[,7:9])
W0lin   <- as.matrix(Ganymede[,10:12])
M0lin   <- as.matrix(Ganymede[,13:15])
W0quadr <- as.matrix(Ganymede[,16:18])
M0quadr <- as.matrix(Ganymede[,19:21])

x0 <- ShowAxialDataSphere(X=X0,V=M0const,vfactor=0.1)
x0 <- ShowAxialDataSphere(X=X0,V=M0lin,vfactor=0.1)
x0 <- ShowAxialDataSphere(X=X0,V=M0quadr,vfactor=0.1)


# Illustrate Bingham distributions ----

na <- 360
th <- 2*pi*(0:na)/na

w <- c(1,1) # 45 deg (=0.78 rad)
w <- c(0,0) # random
w <- c(-2,3)
w <- c(-10,15)

fth <- AxialAngularDensity(w,theta=th)
fth <- AxialAngularHistogramm(w,xymax=2)
# CH: Note that AxialAngularDensity and AxialAngularHistogramm return the same vector fth, 
# and simply plot the data differently

sum(fth[1:na])/na
mean(fth[1:na] >= 1)

kb <- RAlpha(w)
kappa <- kb[1]
sqrt(sum(w^2))
kappa
beta <- kb[2]
beta
w - kappa*c(cos(beta),sin(beta))
r <- Moments1_vMF(kappa,d=2)$M1
r

plot(0.5*(1 + r*cos(2*(th - beta)))*cos(th),
	 0.5*(1 + r*cos(2*(th - beta)))*sin(th),type='l',lwd=2,
	 xlim=c(-1,1),ylim=c(-1,1))
plot(0.5*(1 + cos(2*(th - beta)))*cos(th),
	 0.5*(1 + cos(2*(th - beta)))*sin(th),type='l',lwd=2,
	 xlim=c(-1,1),ylim=c(-1,1))
plot(0.5*(1 + 0.7*cos(2*(th-beta)))*cos(th),
	 0.5*(1 + 0.7*cos(2*(th-beta)))*sin(th),type='l',lwd=2,
	 xlim=c(-1,1),ylim=c(-1,1))


# Illustrate WassersteinBingham ----

# Compare two Bingham distributions:
w0 <- rep(0,2)
r1 <- 3
w1 <- c(r1,0)
r2 <- 2
alpha <- pi/8
w2 <- r2*c(cos(alpha),sin(alpha))
W01 <- WassersteinBingham(w1=w0,w2=w1)
W02 <- WassersteinBingham(w1=w0,w2=w2)
W12 <- WassersteinBingham(w1=w1,w2=w2)
c(W01,W02,W12)
alpha

# Compare a Bingham with another distribution Q:

w <- c(1,2)
beta <- RAlpha(w)[2]
k <- 20
ell=10
q <- rep(dbeta((1:k)/k,ell*beta/pi,ell*(1 - beta/pi)),times=2)
# q[i] is proportional to Q((2*pi*(i-1)/m,2*pi*i/m])
WassersteinBinghamEmpirical(w=w,q=q)
# random: w=c(0,0)
# CH test for main_bingham.R: WassersteinBinghamEmpirical(w=axWs[i,], q=axdens2)



# Illustrate GreedySubset ----

n1 <- 500
X1 <- matrix(rnorm(n1*3),n1,3)
X1 <- X1 / sqrt(rowSums(X1^2))
x0 <- ShowAxialDataSphere(X=X1)

n2 <- 2000
X2 <- matrix(rnorm(n2*3),n2,3)
X2 <- X2 / sqrt(rowSums(X2^2))
X2 <- GreedySubset(X=X2,size=n1)
x0 <- ShowAxialDataSphere(X=X2)

n3 <- 20000
X3 <- matrix(rnorm(n3*3),n3,3)
X3 <- X3 / sqrt(rowSums(X3^2))
X3 <- GreedySubset(X=X3,size=n1)
x0 <- ShowAxialDataSphere(X=X3)


# Test LocalLinQuadrGLM_vMF ----
