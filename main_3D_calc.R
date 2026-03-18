# Caroline Haslebacher, 2025-05-22
# this script smoothes observations with a GLM and fits a vMF distribution for clustered data (handles clusters individually)

#%%
# on CHansle-wks: 
# setwd('D:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# setwd('C:/Users/chaslebacher/Desktop/R_code/Europa_Bingham')
# on algol:
# setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# on laptop:
# setwd('W:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# DELL tower
setwd('E:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')

source("Sunazimuth_IsotonicH.R")

## load libraries
# library(raster)

source('AxialDataSphere.R')
source('GLM_3D.R')
source("Directional.R")

library(foreach)
library(doParallel)
# install with CRAN: install.packages(c("foreach", "doParallel"))

# for test set:
# orpath <- "data/tests/auxiliary_files/for_vMF/"

## specify the data path
# Manual segmentations only (12 files for azimuth study draft)
orpath <- "data/manualsegs/"

# or LM1_1_Regmaps
# orpath <- "for_azimuth_study_draft/LM1_1_Regmaps/"
# orpath <- "data/LM1_1_SSI/"
# LM1_1_Regmpas and LM1_1_SSI combined:
# orpath <- "for_azimuth_study_draft/LM1_1_all/"
# feasibility study
# orpath <- "data/LM1_1_51Peg/"

# or tests:
# orpath <- "for_azimuth_study_draft/test/"

# code adapted from AxialDataSphere_Demo.R
# read the data -
filenames <- Sys.glob(paste(orpath, "*.csv", sep=''))

#%%

# choose either corrected or uncorrected

# sun corrected:
savepath = paste(orpath, "output/", sep='')
suncorrection=TRUE
suncorrection_init <- TRUE
dir.create(savepath, recursive=TRUE, showWarnings=FALSE)

# ## no correction, to compare and interpret results
# savepath = paste(orpath, "output/uncorrected/", sep='')
# suncorrection<-FALSE
# suncorrection_init<-FALSE
# dir.create(savepath, recursive=TRUE, showWarnings=FALSE)


#%%

# TODO: change indexing to start from 0 (and not filenames[4:length(filenames)])
for (filen in filenames) {
	print(filen) # e.g. "for_azimuth_study_draft/manualsegs/ULDR_parents_for_vMF.csv"
	# debug: filen <- "data/manualsegs/RC_children_for_vMF_spherical.csv"
	
	ds <- read.csv(filen,
					header=TRUE,sep=',')
	# construct name to save results
	filesplits <- strsplit(filen, split = '/') # filestem is a nested list. 
	# > filesplits
	# [[1]]
	# [1] "for_azimuth_study_draft"  "manualsegs"               "ULDR_parents_for_vMF.csv"
	# we want to get the last entry
	filestem <- tail(filesplits[[1]], n=1)
	# and then get rid of .csv 
	filestem <- strsplit(filestem, split = '.csv')[[1]] # I can simply strip .csv away! and unpack list
	result_name <- paste(filestem, '_result.csv', sep='')
	rds_result_name <- paste(filestem, '_result.rds', sep='')

	if(grepl('Band', filestem)){
		# setting suncorrection to false for bands
		suncorrection=FALSE
		print('skipping suncorrection for ')
		print(filestem)
	}

	# filter observations by name:
	# observations we keep are
	# 17ESREGMAP02
	# E6ESDRKLIN01
	# 25ESDARKBP01
	# 17ESREGMAP03
	# 17ESREGMAP01EXCERPT2
	# 17ESREGMAP01EXCERPT1
	# 11ESREGMAP01EXCERPT1
	datafiltering <- 	(grepl('17ESREGMAP02', ds$name) | 
						grepl('E6ESDRKLIN01', ds$name) | 
						grepl('25ESDARKBP01', ds$name) | 
						grepl('17ESREGMAP03', ds$name) | 
						grepl('17ESREGMAP01EXCERPT2', ds$name) | 
						grepl('17ESREGMAP01EXCERPT1', ds$name) | 
						grepl('11ESREGMAP01EXCERPT1', ds$name))
	ds <- ds[datafiltering, ]

	# ready the data
	X <- as.matrix(ds[,2:4])
	V <- as.matrix(ds[,8:10])

	# ready the grid:
	C <- ds$cluster # note: clusters are grouped by same emission angle in Generate_input_for_vMF.py.
	# str(C)
	# print(table(C))

	# for example:
	# > table(C)
	# C
	# 0   1   2   3   4   5   6   7   8   9  10  12  14  15  16  17  18
	# 14   4   7   3  49   7   6 599  24  10 246   9  30  32  25 133  14

	# subsetsize <- 200
	# TODO: test a fraction of the number of observations. 20% for example.
	n0 <- sum(as.integer(0.3*table(C)))
	print(as.integer(0.3*table(C)))
	# n0 <- sum(pmin(table(C), subsetsize))# takes whichever is lower, the number of observations in the cluster, or the number 200
	print(paste('n0', n0))

	# CH 2026-01-26: let's make N adaptable. For small clusters, N=200 is too large
	# we calculate N from the number of observations (table C)
	# therefore, we do this inside the for loop further down
	# old code:
	# N <- 200 # we can loop through different N. But for now, let's keep it at 200
	# let's try N=50 too
	a0 <- 0
	set.seed(3032024)

	### initialize data
	X0 <- matrix(0,n0,3)
	V0 <- matrix(0,n0,3)
	W02const <- matrix(0,n0,2) # tmp[1]*c(cos(tmp[2]/2),sin(tmp[2]/2)). The Bingham parameters
	W0const <- matrix(0,n0,3)
	M0const <- matrix(0,n0,3)
	W02lin <- matrix(0,n0,2) # Bingham
	W0lin   <- matrix(0,n0,3)
	M0lin   <- matrix(0,n0,3)
	W02quadr <- matrix(0,n0,2) # Bingham
	W0quadr <- matrix(0,n0,3)
	M0quadr <- matrix(0,n0,3)
	cluster <- matrix(0,n0) # we do not take all values, remember!
	sun <- matrix(0,n0) # keep the subsolar sun azimuth for the bias
	alpha <- matrix(0,n0) # keep alpha as well
	paramN <- matrix(0, n0) # store the variable N parameter (same for one cluster)

	################## customweights calculation
	if(suncorrection){
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

		# get weights by getting the nearest data value from dataframe in R
		# with df[which.min(angle_diff[i],]
		# retrieve the weights easily (pass path to dc):
		# gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3)
		# shift values slightly to match observed curve better:
		if (grepl('DR', filestem)) {
			print('correcting for double ridges')
			# note: grepl is case sensitive, so all good
			csvp <- './isotonic_DTMdatasets/DR_children_for_vMF_spherical.csv'
		}else if (grepl('RC', filestem)) {
			print('correcting for ridge complexes')
			csvp <- './isotonic_DTMdatasets/RC_children_for_vMF_spherical.csv'
		}else if (grepl('UL', filestem)) {
			print('correcting for undifferentiated lineae')
			csvp <- './isotonic_DTMdatasets/UL_children_for_vMF_spherical.csv'
		}else{
			# dff_all correction
			print('correcting for ALL')
			csvp <- './isotonic_DTMdatasets/dff_children_nofilter_for_vMF_spherical.csv' # all
		}

		gsf <- GetSunFit(csvp, M=3, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
		# plot(gsf[,1], gsf[,2])
		# to datafame:
		dfgsf <- data.frame(gsf)
		colnames(dfgsf) <- c('t','ht')
		# retrieve weights:
		customweights <- matrix(0, length(ds[,1]))
		for (i in seq_along(ds[,1])){
			# I am simply finding the correct index with abs(angle_diff[i] - dfgsf$t))
			customweights[i] <- dfgsf[which.min(abs(angle_diff[i] - dfgsf$t)),]$ht
		}
	}else{
		# simply put all customweights to 1
		customweights <- matrix(1, length(ds[,1]))
	}
	print(str(customweights))
	###########

	for (cl in as.numeric(names(table(C)))){ # note: this makes sure we really only take the filled clusters!
		tmp <- which(C==cl)
		# There should not be any empty clusters, but in case:
		if(length(which(C == cl)) <= 1){
			print('empty cluster')
			next
		}
		# get the subset
		tmpsize <- as.integer(0.3*length(ds$name[tmp]))
		if (tmpsize > 0) {
			tmpXV <- GreedySubset(X=X[tmp,],size=tmpsize,V=V[tmp,]) # size=subsetsize
		}else{
			# for example, if there are only 3 points, skip
			# this does not influence a0 and b0, since b0 is 0
			next
		}
		b0 <- dim(tmpXV$X0)[1]
		if (b0 == (tmpsize+1) & (b0 == 2)) {
		   # then, b0 most is 2 and tmpsize is 1, due to a glitch in GreedySubset
		   # we simply only take the first entry and reset b0 to the actual tmpzise
		   b0 <- tmpsize
		   tmpXV$X0 <- tmpXV$X0[1,]
		   tmpXV$V0 <- tmpXV$V0[1,]
		   # ready to continue
		}
		print(paste('Cluster',cl))
		print(paste('size is', tmpsize))
		print(paste('tmpXV dim is', dim(tmpXV$X0)[1]))
		print(paste('a0', a0, ', b0', b0))
		X0[a0 + (1:b0),] <- tmpXV$X0
		V0[a0 + (1:b0),] <- tmpXV$V0
		# # for debugging
		# print(a0)
		# print(a0+b0)
		
		# the below did not work well.
		# # CH 2026-01-26: let's make N adaptable. For small clusters, N=200 is too large
		# # we calculate N from the number of observations (table C)
		# # N should at least be 20, and favorably be a tenth of the number of observations in the cluster
		# if(length(ds$name[tmp]) > 200){
		# 	N <- 100 # 200
		# }else{
		# 	N <- 25 # 50
		# }
		N <- 100 # 50
		# N <- pmax(0.1*length(tmp), 20)
		print(paste('N:', N))
		print(paste('Num observations:', length(ds$name[tmp])))

		tryCatch(
		expr = {
			print('first try')
			regfac <- 0
				res <- SmoothAxialDataSphere(X0=tmpXV$X0, # grid
									X=X[tmp,],V=V[tmp,], # observations
									N=N, regfac=regfac, showsteps=FALSE,
									customweights=customweights[tmp,])
				W02const[a0 + (1:b0),] <- res$W02const
				W0const[a0 + (1:b0),] <- res$W0const
				M0const[a0 + (1:b0),] <- res$M0const
				W02lin[a0 + (1:b0),]   <- res$W02lin
				W0lin[a0 + (1:b0),]   <- res$W0lin
				M0lin[a0 + (1:b0),]   <- res$M0lin
				W02quadr[a0 + (1:b0),] <- res$W02quadr
				W0quadr[a0 + (1:b0),] <- res$W0quadr
				M0quadr[a0 + (1:b0),] <- res$M0quadr
				cluster[a0 + (1:b0),] <- cl
				sun[a0 + (1:b0),] <- ds$sun[tmp][1]
				alpha[a0 + (1:b0),] <- ds$alpha[tmp][1]
				paramN[a0 + (1:b0),] <- N
		},
		error = function(e){ 
			tryCatch(
			expr = {
				print('second try')
				regfac <- 10^(-9)
					res <- SmoothAxialDataSphere(X0=tmpXV$X0,
									X=X[tmp,],V=V[tmp,],
									N=N, regfac=regfac, showsteps=FALSE,
									customweights=customweights[tmp,])
					W02const[a0 + (1:b0),] <- res$W02const
					W0const[a0 + (1:b0),] <- res$W0const
					M0const[a0 + (1:b0),] <- res$M0const
					W02lin[a0 + (1:b0),]   <- res$W02lin
					W0lin[a0 + (1:b0),]   <- res$W0lin
					M0lin[a0 + (1:b0),]   <- res$M0lin
					W02quadr[a0 + (1:b0),] <- res$W02quadr
					W0quadr[a0 + (1:b0),] <- res$W0quadr
					M0quadr[a0 + (1:b0),] <- res$M0quadr
					cluster[a0 + (1:b0),] <- cl
					sun[a0 + (1:b0),] <- ds$sun[tmp][1]
					alpha[a0 + (1:b0),] <- ds$alpha[tmp][1]
					paramN[a0 + (1:b0),] <- N
			},
			error = function(e){
				tryCatch(
				expr = {
					print('third try')
					regfac <- 10^(-5)
						res <- SmoothAxialDataSphere(X0=tmpXV$X0,
									X=X[tmp,],V=V[tmp,],
									N=N, regfac=regfac, showsteps=FALSE,
									customweights=customweights[tmp,])
						W02const[a0 + (1:b0),] <- res$W02const
						W0const[a0 + (1:b0),] <- res$W0const
						M0const[a0 + (1:b0),] <- res$M0const
						W02lin[a0 + (1:b0),]   <- res$W02lin
						W0lin[a0 + (1:b0),]   <- res$W0lin
						M0lin[a0 + (1:b0),]   <- res$M0lin
						W02quadr[a0 + (1:b0),] <- res$W02quadr
						W0quadr[a0 + (1:b0),] <- res$W0quadr
						M0quadr[a0 + (1:b0),] <- res$M0quadr
						cluster[a0 + (1:b0),] <- cl
						sun[a0 + (1:b0),] <- ds$sun[tmp][1]
						alpha[a0 + (1:b0),] <- ds$alpha[tmp][1]
						paramN[a0 + (1:b0),] <- N
				},
				error = function(e){
					print('fourth and last try')
					regfac <- 10^(-1)
					res <- SmoothAxialDataSphere(X0=tmpXV$X0,
									X=X[tmp,],V=V[tmp,],
									N=N, regfac=regfac, showsteps=FALSE,
									customweights=customweights[tmp,])
					W02const[a0 + (1:b0),] <- res$W02const
					W0const[a0 + (1:b0),] <- res$W0const
					M0const[a0 + (1:b0),] <- res$M0const
					W02lin[a0 + (1:b0),]   <- res$W02lin
					W0lin[a0 + (1:b0),]   <- res$W0lin
					M0lin[a0 + (1:b0),]   <- res$M0lin
					W02quadr[a0 + (1:b0),] <- res$W02quadr
					W0quadr[a0 + (1:b0),] <- res$W0quadr
					M0quadr[a0 + (1:b0),] <- res$M0quadr
					cluster[a0 + (1:b0),] <- cl
					sun[a0 + (1:b0),] <- ds$sun[tmp][1]
					alpha[a0 + (1:b0),] <- ds$alpha[tmp][1]
					paramN[a0 + (1:b0),] <- N
				}
				)
			}
			)
		}
		)

		a0 <- a0 + b0
	}


	# TODO: save data
	Europa <- data.frame('X0'=X0,'V0'=V0,
						'W0const'=W0const,'M0const'=M0const,
						'W0lin'=W0lin,'M0lin'=M0lin,
						'W0quadr'=W0quadr,'M0quadr'=M0quadr,
						'cluster'=cluster, 'sun'=sun, 'alpha'=alpha, 'paramN'=paramN, 
						'W02const'=W02const, 'W02lin'=W02lin, 'W02quadr'=W02quadr) # NOTE: we put the Bingham params last, so that we can keep the retrieval as is.
	# delete rows that are zero
	Europa_filtered <- Europa[Europa$paramN != 0, ]
	saveRDS(Europa_filtered,file=paste(savepath, rds_result_name, sep=''))

	if (suncorrection_init) {
		# set it back to TRUE
		suncorrection=TRUE
	}
	# else, we leave it at false (for totally uncorrected)

}



#%%

# # check data:

# # Europa <- readRDS(file="data/manualsegs/output/dff_all_nofilter_for_vMF_spherical_result.rds") # RCB_children_for_vMF_spherical_result.rds")
# # Europa <- readRDS(file="data/manualsegs/output/ULDR_parents_for_vMF_spherical_result.rds")
# # Europa <- readRDS(file="data/manualsegs/output/ULDR_children_for_vMF_spherical_result.rds")
# # Europa <- readRDS(file="data/manualsegs/output/dff_all_nofilter_for_vMF_spherical_result.rds")
# # Europa <- readRDS(file="data/manualsegs/output/dff_children_nofilter_for_vMF_spherical_result.rds")
# Europa <- readRDS(file="data/manualsegs/output/RCB_children_young1_for_vMF_spherical_result.rds")
# # no sun correction:
# Europa <- readRDS(file="data/manualsegs/output_nosun/dff_children_nofilter_for_vMF_spherical_result.rds")

# X0 <- as.matrix(Europa[,1:3])
# V0 <- as.matrix(Europa[,4:6])
# W0const <- as.matrix(Europa[,7:9])
# M0const <- as.matrix(Europa[,10:12])
# W0lin   <- as.matrix(Europa[,13:15])
# M0lin   <- as.matrix(Europa[,16:18])
# W0quadr <- as.matrix(Europa[,19:21])
# M0quadr <- as.matrix(Europa[,22:24])

# # quick bugfix:
# # M0quadr <- as.matrix(Europa[,28:30])

# x0 <- c(1,0,0)
# x0 <- c(0.1930043, -0.9726410, -0.1293014)
# # Local constant fits:
# windows()
# x0 <- ShowAxialDataSphere(X=X0,V=M0const,x0=x0,
# 						  Vcolours = c('blue','gray'))
# # Local linear fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0lin,x0=x0,
# 						  Vcolours = c('blue','gray'))
# # Local quadratic fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0quadr,x0=x0,
# 						  Vcolours = c('blue','gray'))

# #%% auxiliary data to test:
# # check data:
# # initial data:

# # trial1
# # ds <- read.csv("data/tests/auxiliary_files/trial1_initial/for_vMF/dff_all_nofilter_for_vMF_spherical.csv", # _extended
# # 			   header=TRUE,sep=',')
# # str(ds)

# # X <- as.matrix(ds[,2:4])
# # str(X)
# # V <- as.matrix(ds[,8:10])
# # str(V)
# # range(rowSums(X^2))
# # range(rowSums(V^2))
# # range(rowSums(X*V))
# # x0 <- c(1,0,0)
# # x0 <- ShowAxialDataSphere(X=X,V=V, x0=x0,
# #						  Vcolours = c('blue','gray'))


# ds <- read.csv("data/tests/auxiliary_files/for_vMF/dff_all_nofilter_for_vMF_spherical.csv", # _extended
# 			   header=TRUE,sep=',')
# ds <- read.csv("data/tests/auxiliary_files/for_vMF/RCB_parents_for_vMF_spherical.csv", header=TRUE,sep=',')
# # corrected manualsegs:
# ds <- read.csv("data/tests/auxiliary_files/for_vMF/RCB_parents_for_vMF_spherical.csv", header=TRUE,sep=',')
# str(ds)

# X <- as.matrix(ds[,2:4])
# str(X)
# V <- as.matrix(ds[,8:10])
# str(V)
# range(rowSums(X^2))
# range(rowSums(V^2))
# range(rowSums(X*V))
# x0 <- c(1,0,0)
# x0 <- ShowAxialDataSphere(X=X,V=V,x0=x0,
# 						  Vcolours = c('blue','gray'))


# Europa <- readRDS(file="data/tests/auxiliary_files/for_vMF/output/dff_all_nofilter_for_vMF_spherical_result.rds")
# X0 <- as.matrix(Europa[,1:3])
# V0 <- as.matrix(Europa[,4:6])
# W0const <- as.matrix(Europa[,7:9])
# M0const <- as.matrix(Europa[,10:12])
# W0lin   <- as.matrix(Europa[,13:15])
# M0lin   <- as.matrix(Europa[,16:18])
# W0quadr <- as.matrix(Europa[,19:21])
# M0quadr <- as.matrix(Europa[,22:24])
# x0 <- c(1,0,0)
# # show initial data:
# x0 <- ShowAxialDataSphere(X=X0,V=V0,x0=x0,
# 						  Vcolours = c('blue','gray'))

# # Local constant fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0const,x0=x0,
# 						  Vcolours = c('blue','gray'))
# # Local linear fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0lin,x0=x0,
# 						  Vcolours = c('blue','gray'))
# # Local quadratic fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0quadr,x0=x0,
# 						  Vcolours = c('blue','gray'))

# #%%

# ds <- read.csv("ULDR_parents_for_vMF_spherical.csv",
# 			   header=TRUE,sep=',')
# # str(ds)

# X <- as.matrix(ds[,2:4])
# # str(X)
# V <- as.matrix(ds[,8:10])
# # str(V)

# C <- ds$cluster
# str(C)
# table(C)

# subsetsize <- 200
# n0 <- sum(pmin(table(C), subsetsize))# takes whichever is lower, the number of observations in the cluster, or the number 200
# # pmin(table(C), subsetsize) returns
# #   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
# #  56  20  29 108 181 116 118 200 200  74 200  58  45  64 168  50 142 200 155
# # n0 is the sum
# n0

# # Illustrate local weights only for cluster C==7:
# x0 <- ShowWeightedAxialDataSphere(X[C==7,],V[C==7,],N=100, # 200 provides a much better fit than 400. Maybe we can go even lower
# 								  vfactor=0.05)
# x0
# x0 <- ShowWeightedAxialDataSphere(X[C==7,],V[C==7,],N=200,
# 								  vfactor=0.05,x0=x0)
# points(0,0,pch=19,col='red')

# # or for all:
# x0 <- ShowWeightedAxialDataSphere(X,V,N=50, # 200 provides a much better fit than 400. Maybe we can go even lower
# 								  vfactor=0.05)

# ### initialize data
# X0 <- matrix(0,n0,3)
# V0 <- matrix(0,n0,3)
# W0const <- matrix(0,n0,3)
# M0const <- matrix(0,n0,3)
# W0lin   <- matrix(0,n0,3)
# M0lin   <- matrix(0,n0,3)
# W0quadr <- matrix(0,n0,3)
# M0quadr <- matrix(0,n0,3)

# # TODO: parallelize. use code from main_azimuth_calc.R
# N <- 200 # we can loop through different N. Or maybe pick one?
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

# # TODO: save data
# Europa <- data.frame('X0'=X0,'V0'=V0,
# 					'W0const'=W0const,'M0const'=M0const,
# 						'W0lin'=W0lin,'M0lin'=M0lin,
# 						'W0quadr'=W0quadr,'M0quadr'=M0quadr)
# saveRDS(Europa,file='Europa400.rds')


# Europa <- readRDS(file="Europa200.rds")
# X0 <- as.matrix(Europa[,1:3])
# V0 <- as.matrix(Europa[,4:6])
# W0const <- as.matrix(Europa[,7:9])
# M0const <- as.matrix(Europa[,10:12])
# W0lin   <- as.matrix(Europa[,13:15])
# M0lin   <- as.matrix(Europa[,16:18])
# W0quadr <- as.matrix(Europa[,19:21])
# M0quadr <- as.matrix(Europa[,22:24])
# x0 <- c(1,0,0)
# # show initial data:
# x0 <- ShowAxialDataSphere(X=X0,V=V0,x0=x0,
# 						  Vcolours = c('blue','gray'))

# # Local constant fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0const,x0=x0,
# 						  Vcolours = c('blue','gray'))
# # Local linear fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0lin,x0=x0,
# 						  Vcolours = c('blue','gray'))
# # Local quadratic fits:
# x0 <- ShowAxialDataSphere(X=X0,V=M0quadr,x0=x0,
# 						  Vcolours = c('blue','gray'))

# # in a separate script:
# # TODO: for the subgrid, calculate stresses for comparison
# # TODO: bias with sun azi


