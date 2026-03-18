
#### Caroline Haslebacher
# 2026-01-12
# script to plot the histograms of the angle difference (between lineament azimuth and sun azimuth)
# for different subsets (lineament types, clusters)
# 
#%%
# algol:
setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')


#%%

get_sunazi <- function(ds){
    sunazi <- pi/180*(ds$sun)
    # note that sunazi element (-3pi/2, -2pi) because of the transformation
    sunazi[(sunazi >= -pi) & (sunazi < 0)] <- sunazi[(sunazi >= -pi) & (sunazi < 0)] + pi
    sunazi[sunazi < -pi] <- 2*pi + sunazi[sunazi < -pi]
    
    return(sunazi)
}

get_alphs <- function(ds) {
    # make alphas in between 0 and pi (not /pi/2, pi/2)
    alphs <- (pi/2 - ds$alpha) 
    alphs[alphs<0] <- pi + alphs[alphs<0] # e.g. -60deg gets transformed to 180 + (-60deg) = 120deg

    return(alphs)
}

get_anglediff <- function(sunazi, alphs) {
        # angle difference:
    angle_diff <- abs(sunazi - alphs)
    angle_diff[angle_diff>pi/2] <- pi - angle_diff[angle_diff>pi/2]

    return(angle_diff)
}

anglediff_hist <- function(ds, vis_savepath, savename) {
    X <- as.matrix(ds[,2:4])
    V <- as.matrix(ds[,8:10])

    sunazi <- get_sunazi(ds)
    alphs <- get_alphs(ds)
    angle_diff <- get_anglediff(sunazi, alphs)

    # produce and save the hist
    png(filename=paste(vis_savepath, savename, '_anglehist.png', sep=''), width=1500, height=2000, pointsize=30)
    pibins <- seq(0, pi/2, length.out = 16)
    adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias", savename, '\nObservations:', length(angle_diff), sep=' '), xlab="Difference of Sun Azimuth and Lineament Azimuth [rad]", ylab="Probability Density")
    # lines(density(angle_diff, to=pi/2), col='red')
    dev.off()

    # # clusters:
    for(clustername in c('17ESREGMAP02', 'E6ESDRKLIN01', '25ESDARKBP01', '17ESREGMAP03', '17ESREGMAP01EXCERPT2', '17ESREGMAP01EXCERPT1', '11ESREGMAP01EXCERPT1')){
        # print(clustername)
        datafiltering <- (grepl(clustername, ds$name))
        ds_sub <- ds[datafiltering, ]

        sunazi <- get_sunazi(ds_sub)
        alphs <- get_alphs(ds_sub)
        angle_diff <- get_anglediff(sunazi, alphs)

        # produce and save the hist
        png(filename=paste(vis_savepath, savename, '_', clustername, '_anglehist.png', sep=''), width=1500, height=2000, pointsize=30)
        pibins <- seq(0, pi/2, length.out = 16)
        adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias", savename, clustername, '\nObservations:', length(angle_diff), sep=' '), xlab="Difference of Sun Azimuth and Lineament Azimuth [rad]", ylab="Probability Density")
        # lines(density(angle_diff, to=pi/2), col='red')
        dev.off()

    }

    return(adiffh)
}

#%%
base_savepath <- './results/angle_hists/'
dir.create(base_savepath, recursive=TRUE, showWarnings=FALSE)

#%%
orpath <- "data/manualsegs/"
filenames <- Sys.glob(paste(orpath, "*csv", sep=''))

for (fi in seq_along(filenames)) {
    filen <- filenames[fi]
    print(filen) # e.g. "data/manualsegs/output/dff_children_nofilter_for_vMF_spherical_result.rds"
    # construct name to save results
    filesplits <- strsplit(filen, split = '/') # filestem is a nested list. 
    # > filesplits
    # [[1]]
    # [1] "for_azimuth_study_draft"  "manualsegs"               "ULDR_parents_for_vMF.csv"
    # we want to get the last entry
    filestem <- tail(filesplits[[1]], n=1)
    # and then get rid of .csv 
    srcstem <- strsplit(filestem, split = '.csv')[[1]] # I can simply strip .csv away! and unpack list
    savename <- strsplit(srcstem, split = '_for_vMF_spherical')[[1]] # strip away the full '_for_vMF_spherical'

    # MAIN
    ds <- read.csv(filen,
 			   header=TRUE,sep=',')
    adiffh_1 <- anglediff_hist(ds, base_savepath, savename)

    # 

}


#%%


# # only one cluster (or negation)
# cl <- 7
# tmp <- which(C==cl)
# sunazi <- pi/180*(ds$sun[tmp])

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
adiffh_1 <- anglediff_hist(ds, paste(base_savepath, 'dff_all_anglehist.png', sep=''))


# bands:
# parents:
ds <- read.csv("./data/manualsegs/Band_parents_for_vMF_spherical.csv",
 			   header=TRUE,sep=',')
str(ds)
adiffh_1 <- anglediff_hist(ds, paste(base_savepath, 'Band_parents_anglehist.png', sep=''))
# children:
ds <- read.csv("./data/manualsegs/Band_parents_for_vMF_spherical.csv",
 			   header=TRUE,sep=',')
str(ds)
adiffh_1 <- anglediff_hist(ds)


#%% or run:


#%%

#%%
# show sunazimuth hist:
hist(sunazi)
hist(alphs)

#%% actually save it:
# visualisation of the calculated 5623 angle differences (histogram) and the fitting curve I got from the DEM study
# note how well they agree.
vis_savepath <- './results/'


#%% lineament types:








