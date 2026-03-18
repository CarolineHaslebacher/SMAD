
#### Caroline Haslebacher
# 2025-12-03
# 
#%%
# algol:
# setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# DELL tower
setwd('E:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')

source("Sunazimuth_IsotonicH.R")

# load libraries
library(raster)

source('AxialDataSphere.R')
source('GLM_3D.R')
source('Directional.R')




#%%
# We can use AxialAngularDensity(w, theta) and retrieve weights with fth[lineament_azimuth_math_crs]
# for this, let's load some data:
# ds <- read.csv("ULDR_children_for_vMF_spherical.csv",
# 			   header=TRUE,sep=',')
# children
# ds <- read.csv("dff_children_nofilter_for_vMF_spherical_phocube.csv",
# 			   header=TRUE,sep=',')
# PARENTS
ds <- read.csv("./data/manualsegs/dff_all_nofilter_for_vMF_spherical.csv",
 			   header=TRUE,sep=',')
str(ds)

X <- as.matrix(ds[,2:4])
V <- as.matrix(ds[,8:10])


sunazi <- pi/180*(ds$sun)

#
# NEW: use splines weights 1/h(t):
# first, we need to calculate the angle between the sun direction and the azimuth
# sun direction was transformed to mathematical crs in Generate_input_for_vMF.py,
# but alpha was not (therefore, you see pi/2-ds$alpha below)
# symmetry allows to make all sunazimuth values to be inside (0, pi)

# note that sunazi element (-3pi/2, -2pi) because of the transformation
sunazi[(sunazi >= -pi) & (sunazi < 0)] <- sunazi[(sunazi >= -pi) & (sunazi < 0)] + pi
sunazi[sunazi < -pi] <- 2*pi + sunazi[sunazi < -pi]
# make alphas in between 0 and pi (not /pi/2, pi/2)
alphs <- (pi/2 - ds$alpha) 
alphs[alphs<0] <- pi + alphs[alphs<0] # e.g. -60deg gets transformed to 180 + (-60deg) = 120deg

# angle difference:
angle_diff <- abs(sunazi - alphs)
angle_diff[angle_diff>pi/2] <- pi - angle_diff[angle_diff>pi/2]

# get weights by getting the nearest data value from dataframe in R
# with df[which.min(angle_diff[i],]
# retrieve the weights easily (pass path to dc):
# gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3)
# shift a little:
gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
# plot(gsf[,1], gsf[,2])
# to datafame:
dfgsf <- data.frame(gsf)
colnames(dfgsf) <- c('t','ht')
# retrieve weights:
customweights <- matrix(0, length(ds[,1]))
for (i in seq_along(ds[,1])){
    customweights[i] <- dfgsf[which.min(abs(angle_diff[i] - dfgsf$t)),]$ht
}

# look at it:
pibins <- seq(0, pi/2, length.out = 16)
adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias Correction"), xlab="Difference of Sun Azimuth and Lineament Azimuth [rad]", ylab="Probability Density")
lines(dfgsf$t, 1/dfgsf$ht)
# lines(density(angle_diff, to=pi/2), col='red')

# show sunazimuth hist:
hist(sunazi)
hist(alphs)

#%% actually save it:
# visualisation of the calculated 5623 angle differences (histogram) and the fitting curve I got from the DEM study
# note how well they agree.
vis_savepath <- './results/'
png(filename=paste(vis_savepath, 'splines_DEM_WITHdffall_obs_3995PARENTS.png', sep=''), width=1500, height=2000, pointsize=30)
pibins <- seq(0, pi/2, length.out = 16)
# Aik is Difference of Sun Azimuth and Lineament Azimuth [rad]
adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias Correction"), xlab="Aik", ylab="Probability Density")
lines(dfgsf$t, 1/dfgsf$ht, lwd=2, col='red')
# lines(density(angle_diff, to=pi/2), col='red')
dev.off()

vis_savepath <- './results/'
png(filename=paste(vis_savepath, 'splines_DEM_WITHdffall_obs_3995PARENTS_histonly.png', sep=''), width=1500, height=2000, pointsize=30)
pibins <- seq(0, pi/2, length.out = 16)
# Aik is Difference of Sun Azimuth and Lineament Azimuth [rad]
adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias Correction"), xlab="Aik", ylab="Probability Density")
# lines(density(angle_diff, to=pi/2), col='red')
dev.off()



#%%
# and show extracted customweights:
plot(angle_diff, customweights, ylim=c(0,7))
lines(dfgsf$t, 1/dfgsf$ht)
lines(dfgsf$t, dfgsf$ht)


#%% lineament categories (parameter id_int)
# id_int = 1 are bands, but there are none in the DTM region
# id_int = 2 are double ridges
# id_int = 3 are ridge complexes
# id_int = 4 are undifferentiated lineae

# # for debugging, choose one, comment out the others:
# # all:
# savestr <- 'All lineaments'
# # double ridges
# savestr <- 'Double Ridges'
# # ridge complexes
# savestr <- 'Ridge Complexes'
# # undifferentiated lineae
# savestr <- 'Undifferentiated Lineae'

cat_names <- c('All lineaments', 'Double Ridges', 'Ridge Complexes', 'Undifferentiated Lineae')
# select rank of splines
Msel <- 3

for(savestr in cat_names){
	if(savestr == 'All lineaments'){
	csvp <- './isotonic_DTMdatasets/dff_children_nofilter_for_vMF_spherical.csv' # all
	ds <- read.csv("./data/manualsegs/dff_all_nofilter_for_vMF_spherical.csv",
				header=TRUE,sep=',')
	}else if (savestr == 'Double Ridges') {
		csvp <- './isotonic_DTMdatasets/DR_children_for_vMF_spherical.csv'
		ds <- read.csv("./data/manualsegs/DR_children_for_vMF_spherical.csv",
					header=TRUE,sep=',')
	}else if (savestr == 'Ridge Complexes') {
		csvp <- './isotonic_DTMdatasets/RC_children_for_vMF_spherical.csv'
		ds <- read.csv("./data/manualsegs/RC_children_for_vMF_spherical.csv",
					header=TRUE,sep=',')
	}else if (savestr == 'Undifferentiated Lineae') {
		csvp <- './isotonic_DTMdatasets/UL_children_for_vMF_spherical.csv'
		ds <- read.csv("./data/manualsegs/UL_children_for_vMF_spherical.csv",
					header=TRUE,sep=',')
	}

	str(ds)
	X <- as.matrix(ds[,2:4])
	V <- as.matrix(ds[,8:10])
	sunazi <- pi/180*(ds$sun)

	# note that sunazi element (-3pi/2, -2pi) because of the transformation
	sunazi[(sunazi >= -pi) & (sunazi < 0)] <- sunazi[(sunazi >= -pi) & (sunazi < 0)] + pi
	sunazi[sunazi < -pi] <- 2*pi + sunazi[sunazi < -pi]
	# make alphas in between 0 and pi (not /pi/2, pi/2)
	alphs <- (pi/2 - ds$alpha) 
	alphs[alphs<0] <- pi + alphs[alphs<0] # e.g. -60deg gets transformed to 180 + (-60deg) = 120deg

	# angle difference:
	angle_diff <- abs(sunazi - alphs)
	angle_diff[angle_diff>pi/2] <- pi - angle_diff[angle_diff>pi/2]

	gsf <- GetSunFit(csvp, M=Msel, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
	# plot(gsf[,1], gsf[,2])
	# to datafame:
	dfgsf <- data.frame(gsf)
	colnames(dfgsf) <- c('t','ht')
	# retrieve weights:
	customweights <- matrix(0, length(ds[,1]))
	for (i in seq_along(ds[,1])){
		customweights[i] <- dfgsf[which.min(abs(angle_diff[i] - dfgsf$t)),]$ht
	}
	# for the first time, we initialize the data frame
	if (savestr == 'All lineaments') {
	   corr_funcs <- dfgsf
	}else{
		# we only append the ht function:
		corr_funcs <- cbind(corr_funcs, dfgsf$ht)
	}

	######### plot A,J original visibility study data as well
	dc <- read.csv(csvp)
	A <- c(dc$A00, dc$A45, dc$A90, dc$A135)
	J <- c(dc$J00, dc$J45, dc$J90, dc$J135)
	# modify so that Aik element [0, pi/2] instead of [0, pi]
	A[A>pi/2] <- pi - A[A>pi/2]

	# shift if desired:
	Jshifta <- 0.1 # former 0
	Jshiftb <- 0.4 # former 0.33
	Jshiftc <- 0.7 # former 0.67
	# no shfit for 1
	J[J==0] <- Jshifta
	J[J==0.33] <- Jshiftb
	J[J==0.67] <- Jshiftc
	####### ready to be plotted

	adiffh <- hist(angle_diff, freq=TRUE, breaks=pibins)
	# rescale to 0-1
	adiffh$counts <- adiffh$counts/max(adiffh$counts)

	# save it
	vis_savepath <- './results/'
	png(filename=paste(vis_savepath, 'splines_Msel', Msel ,'_DTM_', savestr, '.png', sep=''), width=1500, height=2000, pointsize=40)
	plot(adiffh, main=savestr, xlab="Aik", ylab="Probability Density (scaled to 0-1)", xlim=c(0, pi/2))
	axis(2, at = c(0.2, 0.4, 0.6, 0.8))
	points(A, J)
	pibins <- seq(0, pi/2, length.out = 16)
	lines(dfgsf$t, 1/dfgsf$ht, lwd=2, col='red')
	# lines(density(angle_diff, to=pi/2), col='red')
	dev.off()
}
colnames(corr_funcs) <- c('x', cat_names)
linestyles <- c('solid', 'dashed', "dotted", 'twodash')

png(filename=paste(vis_savepath, 'correction_functions.png', sep=''), width=1500, height=2000, pointsize=40)
# plot all correction functions:
plot(corr_funcs$x, corr_funcs[,2], type='l', lty=linestyles[1], ylim=c(0,max(corr_funcs)), lwd=3,
	main="Sun Bias Correction Functions", ylab="Weights (1/h(Aik))", xlab="Angular Difference (Aik)")
lines(corr_funcs$x, corr_funcs[,3], lty=linestyles[2], lwd=3)
lines(corr_funcs$x, corr_funcs[,4], lty=linestyles[3], lwd=3)
lines(corr_funcs$x, corr_funcs[,5], lty=linestyles[4], lwd=3)
# Add the legend
legend("topright",               # Position keyword
       cat_names, # Vector of labels
       col = 'black',       # Vector of corresponding colors
       lty = linestyles,              # Vector of corresponding line types (1 = solid)
       lwd = c(3,3,3,3)              # Vector of corresponding line widths
       )
dev.off()

#%%
# plot bands without any DTM splines data (there aren't any bands in the DTM)
savestr <- 'Bands'
ds <- read.csv("./data/manualsegs/Band_children_for_vMF_spherical.csv",
				header=TRUE,sep=',')
str(ds)
X <- as.matrix(ds[,2:4])
V <- as.matrix(ds[,8:10])
sunazi <- pi/180*(ds$sun)

# note that sunazi element (-3pi/2, -2pi) because of the transformation
sunazi[(sunazi >= -pi) & (sunazi < 0)] <- sunazi[(sunazi >= -pi) & (sunazi < 0)] + pi
sunazi[sunazi < -pi] <- 2*pi + sunazi[sunazi < -pi]
# make alphas in between 0 and pi (not /pi/2, pi/2)
alphs <- (pi/2 - ds$alpha) 
alphs[alphs<0] <- pi + alphs[alphs<0] # e.g. -60deg gets transformed to 180 + (-60deg) = 120deg

# angle difference:
angle_diff <- abs(sunazi - alphs)
angle_diff[angle_diff>pi/2] <- pi - angle_diff[angle_diff>pi/2]

adiffh <- hist(angle_diff, freq=TRUE, breaks=pibins)
# rescale to 0-1
adiffh$counts <- adiffh$counts/max(adiffh$counts)

vis_savepath <- './results/'
png(filename=paste(vis_savepath, 'histogram_', savestr, '.png', sep=''), width=1500, height=2000, pointsize=40)
plot(adiffh, main=savestr, xlab="Aik", ylab="Probability Density (scaled to 0-1)", xlim=c(0, pi/2))
axis(2, at = c(0.2, 0.4, 0.6, 0.8))
pibins <- seq(0, pi/2, length.out = 16)
dev.off()

#%%


