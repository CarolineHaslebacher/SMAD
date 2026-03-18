
#%%
# algol:
# setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# DELL tower
setwd('E:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
source("Sunazimuth_IsotonicH.R")


#%% CH check data
dc <- read.csv('dff_children_nofilter_for_isotonic_spherical.csv')
# > str(dc)
# 'data.frame':   242 obs. of  24 variables:
#  $ X               : int  0 1 2 3 4 5 6 7 8 9 ...
#  $ X1              : num  -0.565 -0.577 -0.573 -0.577 -0.566 ...
#  $ X2              : num  0.569 0.58 0.584 0.574 0.581 ...
#  $ X3              : num  0.598 0.574 0.575 0.581 0.584 ...
#  $ Y1              : num  -0.693 -0.813 -0.816 -0.695 0.597 ...
#  $ Y2              : num  0.0661 -0.4698 -0.3444 0.0283 0.7782 ...
#  $ Y3              : num  -0.718 -0.343 -0.464 -0.718 -0.196 ...
#  $ V1              : num  -0.594 -0.379 -0.444 -0.584 -0.815 ...
#  $ V2              : num  -0.783 -0.814 -0.811 -0.787 -0.292 ...
#  $ V3              : num  0.183 0.441 0.381 0.197 -0.5 ...
#  $ alpha           : num  1.34 1 1.09 1.33 2.23 ...
#  $ lon             : num  2.35 2.35 2.35 2.36 2.34 ...
#  $ lat             : num  0.641 0.612 0.613 0.62 0.624 ...
#  $ scalarproduct_XY: num  0 0 0 0 0 0 0 0 0 0 ...
#  $ sun             : num  -0.092 -0.092 -0.092 -0.092 -0.092 ...
#  $ cluster         : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ J00             : num  1 0.67 1 0.67 0.33 1 0.33 0.33 1 0.33 ...
#  $ J45             : num  0.33 0 0.33 0.33 1 0.67 0 0.67 0.67 0 ...
#  $ J90             : num  0 0.33 0.67 0.67 0.67 0 0.67 1 0 0.67 ...
#  $ J135            : num  0.67 1 1 1 0.33 0.67 1 1 0.67 1 ...
#  $ A00             : num  1.34 1 1.09 1.33 2.23 ...
#  $ A45             : num  0.555 0.216 0.301 0.541 1.449 ...
#  $ A90             : num  0.231 0.569 0.484 0.245 0.664 ...
#  $ A135            : num  1.016 1.355 1.269 1.03 0.122 ...

# check
plot(dc$A00, dc$J00)
points(dc$A45, dc$J45)
points(dc$A90, dc$J90)
points(dc$A135, dc$J135)

A <- c(dc$A00, dc$A45, dc$A90, dc$A135)
J <- c(dc$J00, dc$J45, dc$J90, dc$J135)
# Check:
plot(A, J)



# modify so that Aik element [0, pi/2] instead of [0, pi]
A[A>pi/2] <- pi - A[A>pi/2]
# check:
plot(A, J)

# shift if desired:
Jshifta <- 0.1 # former 0
Jshiftb <- 0.4 # former 0.33
Jshiftc <- 0.7 # former 0.67
# no shfit for 1
J[J==0] <- Jshifta
J[J==0.33] <- Jshiftb
J[J==0.67] <- Jshiftc
# check:
plot(A, J, ylim=c(0,1))


N <- length(dc$J45)
m <- 4
AJ <- matrix(0,N,2*m) # 2*m, 1 for A, 1 for J
AJ[,seq(1,2*m-1,2)] <- A
AJ[,seq(2,2*m,2)] <- J

plot(AJ)

plot(A,J,xlab='t',ylab='h(t)', ylim=c(0,1))

# Isotonic fit:
TH <- IsotonicH(AJ)
lines(TH[,1],TH[,2],lwd=2,col='black')

# Parametric fit:
XH <- FittedH(AJ,M=6)
lines(XH[,1],XH[,2],lwd=2,col='blue')
XH <- FittedH(AJ,M=4)
lines(XH[,1],XH[,2],lwd=2,col='red')
XH <- FittedH(AJ,M=3)
lines(XH[,1],XH[,2],lwd=2,col='green')
XH <- FittedH(AJ,M=2)
lines(XH[,1],XH[,2],lwd=2,col='orange') # ends at 1, but starts above

xs <- seq(0, pi/2, 0.01)
lines(xs,sin(xs), lwd=2, col='purple')

# final fit:
plot(A,J,xlab='t',ylab='h(t)', ylim=c(0,1))
# Isotonic fit:
TH <- IsotonicH(AJ)
lines(TH[,1],TH[,2],lwd=2,col='black')
# Parametric fit:
XH <- FittedH(AJ,M=3)
lines(XH[,1],XH[,2],lwd=2,col='green')


## show final weighting function:
# 1/h(t)
XH <- FittedH(AJ,M=3)
plot(XH[,1],1/XH[,2],lwd=2,col='green')
points(A,J,xlab='t',ylab='h(t)')
# Isotonic fit:
TH <- IsotonicH(AJ)
lines(TH[,1],TH[,2],lwd=2,col='black')

## finally, I've implemented a function:
# retrieve the weights easily (pass path to dc):
source("Sunazimuth_IsotonicH.R")
gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
plot(gsf[,1], gsf[,2])
plot(gsf[,1], 1/gsf[,2], ylim=c(0,1)) #'original curve'


#%% figures for publication:

png(filename='./results/M3_Splines_fit_visibilitycurve.png', width = 1500, height = 1500, units = "px", pointsize = 32)
# final fit:
# expression(A[i]^k)
plot(A,J,xlab='Angular Difference Aik',ylab='Visibility h(Aik)', ylim=c(0,1), pch = 21, bg = rgb(0, 0.5, 1, 0.2), col = 'black', # Blue color with 50% transparency (alpha = 0.5)
 cex = 2)
# Parametric fit:
XH <- FittedH(AJ,M=3)
lines(XH[,1],XH[,2],lwd=3,col='red')
legend("topleft", legend = c(expression(h(A[i]^k)), expression(J[i]^k)), col = c("red", 'black'), 
        lty = c(1, NA), lwd = c(3, NA), pch = c(NA, 21), pt.cex=2, pt.bg=rgb(0, 0.5, 1, 0.2), cex = 1) # cex is size
dev.off()


# 1/h(t)
png(filename='./results/M3_Splines_fit_sunbiasfunction.png', width = 1500, height = 1500, units = "px", pointsize = 32)
# 1/h(t) final sun bias function
# Parametric fit:
plot(A,J,xlab='Angular Difference Aik',ylab='Weights 1/h(Aik)', ylim=c(0,7), pch = 21, bg = rgb(0, 0.5, 1, 0.2), col = 'black', # Blue color with 50% transparency (alpha = 0.5)
 cex = 1)
XH <- FittedH(AJ,M=3)
points(XH[,1],1/XH[,2], type='l',lwd=3,col='#c49a31')
lines(XH[,1],XH[,2],lwd=3,col='red')
legend("topright", legend = c(expression(1/h(A[i]^k)), expression(h(A[i]^k)), expression(J[i]^k)), col = c('#c49a31', 'red', 'black'), 
        lty = c(1, 1, NA), lwd = c(3, 3, NA), pch = c(NA, NA, 21), pt.cex=1, pt.bg=rgb(0, 0.5, 1, 0.2), cex = 1)
dev.off()


