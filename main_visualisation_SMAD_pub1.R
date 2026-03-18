# 2025-07-09
# Caroline Haslebacher
# this script plots the seven regions selected for closer analysis with SMAD

#%%
# on CHansle-wks: 
# setwd('D:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# mounted algol storage:
# setwd('W:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# setwd('C:/Users/chaslebacher/Desktop/R_code/Europa_Bingham')
# on algol:
# setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# setwd('B:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')

# DELL tower
setwd('E:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')

# load libraries
library(raster)

source('AxialDataSphere.R')
source('GLM_3D.R')
source('Directional.R')

source("Sunazimuth_IsotonicH.R")

library(foreach)
library(doParallel)

set.seed(3032024)



#%%
main <- function(filenames, vis_savepath_base, location=NULL, x0_list=NULL, cluster=TRUE, suncorrection=TRUE, zoomfactor=1, plot=TRUE, rawdata=TRUE, dimension=1500, color='blue', lwd=4){

    # We can either use location and x0_list, for example for the hemispheres,
    # or pass cluster=TRUE to plot each cluster separately at the center

    # for debugging/testing:
    # location=location
    # x0_list=x0_list
    # cluster=FALSE
    # zoomfactor=2**(1/4) (Standard value in AxialDataSphere.R)
    # plot=TRUE

    df <- data.frame()

    # where to save
    # make directory
    dir.create(vis_savepath_base, recursive=TRUE, showWarnings=FALSE)
    # loop through files
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
        srcstem <- strsplit(filestem, split = '.rds')[[1]] # I can simply strip .csv away! and unpack list

        # load real V for the same grid:
        Europa <- readRDS(file=filen)
        # Europa <- readRDS(file="data/manualsegs/output/dff_all_nofilter_for_vMF_spherical.rds")
        X0_Europa <- as.matrix(Europa[,1:3])
        V0 <- as.matrix(Europa[,4:6])
        W0const <- as.matrix(Europa[,7:9])
        M0const <- as.matrix(Europa[,10:12])
        W0lin   <- as.matrix(Europa[,13:15])
        M0lin   <- as.matrix(Europa[,16:18])
        W0quadr <- as.matrix(Europa[,19:21])
        M0quadr <- as.matrix(Europa[,22:24])

        C <- Europa$cluster # note: clusters are grouped by same emission angle in Generate_input_for_vMF.py.
        # table(C)
        # subgrid:
        # we can use the same function to retrieve the indexes: (adapted)
        # we randomly select exactly size=200 grid points on the whole sphere
        size_displaygrid <- 150
        if (length(Europa[,1]) >= size_displaygrid) {
            tmpXV <- GreedySubset(X=X0_Europa,size=size_displaygrid, returnJ=TRUE)
            X0_grid <- tmpXV$X0
            J <- tmpXV$J # these are the indexes
        }else {
            # use full grid (otherwise the GreedySubset can easily break)
            X0_grid <- X0_Europa
            J <- 1:length(Europa[,1])
        }
        # drop J values for (0, 0, 0) X0 values (there are a few cases where X0 is 0 (maybe populated from (0,0,0,0) empty cases from the MAIN loops (could be fixed there as well)))
        new_J <- list()
        for(jitem in J){
            # print(jitem)
            if((X0_Europa[jitem,1] == 0) && (X0_Europa[jitem,1] == 0) && (X0_Europa[jitem,1] == 0)){
                print('culprit found')
            }else{
                new_J <- c(new_J, jitem)
            }
        }
        # assign cleaned J:
        J <- unlist(new_J)

        ### plotting. zoom into each cluster.
        if(cluster){

            # for the weights plot, use full initial data. load it
            filen_csv_orig <- paste("data/manualsegs/", strsplit(srcstem, '_result'),'.csv', sep='')
            ds <- read.csv(filen_csv_orig,
					header=TRUE,sep=',')
            # grepl(substring, string) checks if substring is contained in string
            datafiltering <- (grepl('17ESREGMAP02', ds$name) | # cluster 7
                                grepl('E6ESDRKLIN01', ds$name) |  # cluster 10
                                grepl('25ESDARKBP01', ds$name) | # cluster 8
                                grepl('17ESREGMAP03', ds$name) |  # cluster 17
                                grepl('17ESREGMAP01EXCERPT2', ds$name) |  # cluster 16
                                grepl('17ESREGMAP01EXCERPT1', ds$name) |  # cluster 18
                                grepl('11ESREGMAP01EXCERPT1', ds$name)) # cluster 4
            ds <- ds[datafiltering, ]
            # save number of observations for this filename
            new_file_df <- as.data.frame(table(ds$cluster))
            colnames(new_file_df) <- c("cluster", strsplit(srcstem, split = '_for_vMF_spherical_result')[[1]])
            rownames(new_file_df) <- c("11ESREGMAP01EXCERPT1", # cluster 4  
                                        "17ESREGMAP02", # cluster 7
                                        "25ESDARKBP01", # cluster 8
                                        "E6ESDRKLIN01", # cluster 10
                                        "17ESREGMAP01EXCERPT2", # cluster 16
                                        "17ESREGMAP03", # cluster 17
                                        "17ESREGMAP01EXCERPT1") # cluster 18


            if("cluster" %in% colnames(df)){
                # do not append the 'cluster' column
                df <- cbind(df, new_file_df[2])
            }else{
                # initiate the 'cluster' column and use rbind for initialization of the rows
                # CATCH: this only works if the first file has all rows (no missing observations)!
                df <- rbind(df, new_file_df)
            }

            if(plot){
                X_orig <- as.matrix(ds[,2:4])
                V_orig <- as.matrix(ds[,8:10])
                if(suncorrection){
                    ################## customweights calculation - prepare ShowWeightedAxialDataSphere
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
                        customweights[i] <- dfgsf[which.min(abs(angle_diff[i] - dfgsf$t)),]$ht
                    }
                    ###########
                }else {
                customweights <- matrix(1, length(ds[,1]))
                }
                
                for(i in as.integer(names(table(C)))){

                        print(i)

                        tmp <- which(C==i)

                        if(length(tmp)*Europa[tmp,]$paramN[1] > 0){
                            # note: both length(tmp) and Europa[tmp,]$paramN[1] have to be greater than zero.
                            # I realized that Europa[tmp,]$paramN[1] is 0 if there is an empty cluster (I investigated this).

                            # determine x0 as a selected point in the cluster
                            x0_cluster <- c(mean(X0_Europa[tmp,][,1]), mean(X0_Europa[tmp,][,2]), mean(X0_Europa[tmp,][,3])) # selects approximately the middle point
                            print(x0_cluster)
                            # earlier: # X0_Europa[tmp,][length(tmp)/2,] # random point
                            # > X0_Europa[tmp,][100,]
                            #     X0.1       X0.2       X0.3 
                            # 0.1510836 -0.6963764 -0.7015936


                            # create one subdirectory per location of interest (this makes it much easier to compare different datasets, as the plots in one folder show the same projection)
                            # dir.create(paste(vis_savepath_base, '/', location[i], '/', sep=''), recursive=TRUE, showWarnings=TRUE)

                            # GLM model results
                            # initial GLM output:
                            model_attribute <- 'GLMquadr'
                            # special case for cluster 7 (17ESREGMAP02): increase vfactor
                            if (i == 7) {
                               vfactor <- 0.1
                            }else{
                                vfactor <- 0.05
                                # next # quick fix to try this out
                            }
                            png(filename=paste(vis_savepath_base, '/', srcstem, '_', 'cluster', i, '_', model_attribute, '.png', sep=''), width=dimension, height=dimension)
                            ShowAxialDataSphere(x0=x0_cluster, X=X0_Europa, V=M0quadr, vfactor=vfactor, zoomfactor=zoomfactor,
                                Vcolours=c(color, 'grey'), basemapImage=grey_basemap, lwdVfactor=lwd,
                                use_locator=FALSE)
                            dev.off()

                            if (rawdata) {
                                # initial raw observations on grid:
                                model_attribute <- 'obs'
                                png(filename=paste(vis_savepath_base, '/', srcstem, '_', 'cluster', i, '_', model_attribute, '.png', sep=''), width=dimension, height=dimension)
                                ShowAxialDataSphere(x0=x0_cluster, X=X0_Europa, V=V0, vfactor=0.05, zoomfactor=zoomfactor,
                                    Vcolours=c(color, 'grey'), basemapImage=grey_basemap, lwdVfactor=1,
                                    use_locator=FALSE)
                                dev.off()
                            }

                            # weighted (using original data)
                            model_attribute <- 'Weights'
                            png(filename=paste(vis_savepath_base, '/', srcstem, '_', 'cluster', i, '_', model_attribute, '.png', sep=''), width=dimension, height=dimension)
                            ShowWeightedAxialDataSphere(X_orig, V_orig, x0=x0_cluster, zoomfactor=zoomfactor, N=Europa[tmp,]$paramN[1], # retrieve parameter N here
                                        vfactor=0.05, customweights=customweights, use_locator=FALSE)
                            dev.off()
                        }else{
                            print('have to skip')
                            print(length(tmp))
                        }


                    }
            }
        }else{
                ### plotting
                if(plot){
                    for (i in 1:length(location)){
                        print(i)

                        # create one subdirectory per location of interest (this makes it much easier to compare different datasets, as the plots in one folder show the same projection)
                        # dir.create(paste(vis_savepath_base, '/', location[i], '/', sep=''), recursive=TRUE, showWarnings=TRUE)

                        # GLM model results
                        # initial GLM output:
                        model_attribute <- 'GLMquadr'
                        png(filename=paste(vis_savepath_base, '/', srcstem, '_AND_', model_attribute, '_', location[i], '.png', sep=''), width=dimension, height=dimension)
                        ShowAxialDataSphere(x0=x0_list[,i], X=X0_Europa, V=M0quadr, vfactor=0.05, zoomfactor=zoomfactor,
                            Vcolours=c(color, 'grey'), basemapImage=grey_basemap, lwdVfactor=lwd,
                            use_locator=FALSE)
                        dev.off()

                        # initial raw observations on grid:
                        if (rawdata) {
                            model_attribute <- 'obs'
                            png(filename=paste(vis_savepath_base, '/', srcstem, '_AND_', model_attribute, '_', location[i], '.png', sep=''), width=dimension, height=dimension)
                            ShowAxialDataSphere(x0=x0_list[,i], X=X0_Europa, V=V0, vfactor=0.05, zoomfactor=zoomfactor,
                                Vcolours=c(color, 'grey'), basemapImage=grey_basemap, lwdVfactor=1,
                                use_locator=FALSE)
                            dev.off()
                           
                        }

                    }
                }
        }


    }
    return(df)
}

#%% MAIN
# define centers for visualisation (based on clusters)
# 3 defined x0 to show all 7 regions
location <- c('cl16', 'cl7', 'cl10', 'cl8')
# generate a matrix (I call it 'x0_list') with the locations to display (x0)
x0_cl16 <- c(-0.4031127, 0.2078205, -0.8907541)
x0_cl7<- c(0.13, -0.65, -0.535)
x0_cl10 <- c(-0.3, 0.95, 0.11)
x0_cl8 <- c(0.73, 0.09, 0.68)

x0_list <- cbind(x0_cl16, x0_cl7, x0_cl10, x0_cl8)


# else, we run the code for all possible combinations and produce visualisation plots for both hemispheres
# 

# failsave for foreach
# gc()  # Perform garbage collection to free up memory/connection resources
# options(connection.limit = 500)  # Increase the connection limit if needed


# # 4 hemispheres
# location <- c('leadSubj', 'trailSubj', 'leadAntij', 'trailAntij')
# # generate a matrix (I call it 'x0_list') with the locations to display (x0)
# x0_lead_subj <- c(1, -0.15, 0)
# x0_trail_subj <- c(0.64, 0.77, 0)
# x0_lead_antij <- c(-0.13, -1, 0)
# x0_trail_antij <- c(-0.91, 0.39, 0)
# x0_list <- cbind(x0_lead_subj, x0_trail_subj, x0_lead_antij, x0_trail_antij)
# # moved inside oblqDeg iteration # create one subdirectory per location of interest (this makes it much easier to compare different datasets, as the plots in one folder show the same projection)
# # for (loc in location) {
# #     dir.create(paste(vis_savepath, '/', loc, '/', sep=''), recursive=TRUE, showWarnings=TRUE)
# # }


#%%

# load basemap basemap:
# grey_basemap <- raster("GEOREF_Europa_1200.tif")
grey_basemap <- raster("Basemap_Europa_6LS_3600.tif") # 1200, 2400, or 3600 is available
# Ganymede:
# grey_basemap <- raster("Ganymede_basemap_1200.tif")
# radius <- 1560800 # m
extent(grey_basemap) <- c(-180, 180, -90, 90) # set extent

#%% run all
cluster_based <- TRUE

#%%
# # manual segs, sun corrected, full size
suncorrection <- TRUE
orpath <- "data/manualsegs/output/"
####
zoomfactor <- 1.2
lwd <- 1
vis_savepath_base <- './results/visualisation/preferred_directions/'
filenames <- Sys.glob(paste(orpath, "*.rds", sep='')) # if only used for dff_all, use: filenames <- filenames[9:10]
### main:
main(filenames, vis_savepath_base, location=location, x0_list=x0_list, cluster=FALSE, rawdata=TRUE, dimension=1500, color='blue', zoomfactor=zoomfactor, lwd=lwd)
if(cluster_based){
    # # cluster-based
    vis_savepath_base_clusters <- paste(vis_savepath_base, 'clusters/', sep='')
    main(filenames, vis_savepath_base_clusters, cluster=TRUE, suncorrection = suncorrection, rawdata=TRUE, dimension=1500, color='blue')
}

#%%
# # manual segs, sun corrected, zoomfactor = 2
suncorrection <- TRUE
orpath <- "data/manualsegs/output/"
###
zoomfactor <- 2
lwd<-4
vis_savepath_base <- './results/visualisation/preferred_directions/zoomfactor/'
filenames <- Sys.glob(paste(orpath, "*.rds", sep='')) # if only used for dff_all, use: filenames <- filenames[9:10]
### main:
main(filenames, vis_savepath_base, location=location, x0_list=x0_list, cluster=FALSE, rawdata=TRUE, dimension=1500, color='blue', zoomfactor=zoomfactor, lwd=lwd)
if(cluster_based){
    # clusters with zoomfactor=2.5:
    vis_savepath_base <- './results/visualisation/preferred_directions/zoomfactor/'
    vis_savepath_base_clusters <- paste(vis_savepath_base, 'clusters/', sep='')
    main(filenames, vis_savepath_base_clusters, cluster=TRUE, suncorrection = suncorrection, zoomfactor=2.5, rawdata=FALSE, dimension=1500, color='blue')
}

#%%

# manual segs, uncorrected, full size
suncorrection <- FALSE
orpath <- "data/manualsegs/output/uncorrected/"
zoomfactor <- 1.2
lwd <- 1
vis_savepath_base <- './results/visualisation/preferred_directions/uncorrected/'
filenames <- Sys.glob(paste(orpath, "*.rds", sep='')) # if only used for dff_all, use: filenames <- filenames[9:10]
### main:
main(filenames, vis_savepath_base, location=location, x0_list=x0_list, cluster=FALSE, rawdata=TRUE, dimension=1500, color='blue', zoomfactor=zoomfactor, lwd=lwd)
if(cluster_based){
    # # cluster-based
    vis_savepath_base_clusters <- paste(vis_savepath_base, 'clusters/', sep='')
    main(filenames, vis_savepath_base_clusters, cluster=TRUE, suncorrection = suncorrection, rawdata=TRUE, dimension=1500, color='blue')
}

#%%
# manual segs, uncorrected, zoomfactor=2
suncorrection <- FALSE
orpath <- "data/manualsegs/output/uncorrected/"
zoomfactor <- 2
lwd<-4
vis_savepath_base <- './results/visualisation/preferred_directions/zoomfactor/uncorrected/'
filenames <- Sys.glob(paste(orpath, "*.rds", sep='')) # if only used for dff_all, use: filenames <- filenames[9:10]
### main:
main(filenames, vis_savepath_base, location=location, x0_list=x0_list, cluster=FALSE, rawdata=TRUE, dimension=1500, color='blue', zoomfactor=zoomfactor, lwd=lwd)
if(cluster_based){
    # clusters with zoomfactor=2.5:
    vis_savepath_base <- './results/visualisation/preferred_directions/zoomfactor/'
    vis_savepath_base_clusters <- paste(vis_savepath_base, 'clusters/', sep='')
    main(filenames, vis_savepath_base_clusters, cluster=TRUE, suncorrection = suncorrection, zoomfactor=2.5, rawdata=FALSE, dimension=1500, color='blue')
}

#%%
# print a table of number of observations as well:
# df only:
df <- main(filenames, vis_savepath_base_clusters, cluster=TRUE, plot=FALSE)
# save dataframe to csv
write.csv(df, file = "./results/visualisation/preferred_directions/Number_of_Observations.csv", row.names = TRUE)

#%%
# # for AGU, define filenames manually
# filenames <- c(
#     "data/manualsegs/output/dff_all_nofilter_for_vMF_spherical_result.rds",
#     "data/manualsegs/output/Band_parents_for_vMF_spherical_result.rds",
#     "data/manualsegs/output/RC_parents_for_vMF_spherical_result.rds",
#     "data/manualsegs/output/DR_parents_for_vMF_spherical_result.rds",
#     "data/manualsegs/output/UL_parents_for_vMF_spherical_result.rds"
# )
#%%



#%%
# # 4 hemispheres
# main(filenames, vis_savepath_base, location=location, x0_list=x0_list, cluster=FALSE)


#%%
# developing WhBgh distance


#%%
# source file stem
# choose one:
# srcstem <- "dff_all_nofilter_for_vMF_spherical"
# srcstem <- "RCB_parents_for_vMF_spherical"
# srcstem <- "RCB_children_young2_for_vMF_spherical"
# srcstem <- "RCB_children_old2_for_vMF_spherical"
# srcstem <- "ULDR_parents_for_vMF_spherical"
# srcstem <- "ULDR_children_young2_for_vMF_spherical"
# srcstem <- "ULDR_children_old2_for_vMF_spherical"
# srcstem_list <- cbind(
#     "dff_all_nofilter_for_vMF_spherical",
#     'dff_children_nofilter_for_vMF',
#     'RCB_children_for_vMF',
#     'RCB_children_old1_for_vMF',
#     "RCB_children_old2_for_vMF_spherical",
#     'RCB_children_young1_for_vMF',
#     "RCB_children_young2_for_vMF_spherical",
#     "RCB_parents_for_vMF_spherical",
#     "ULDR_children_for_vMF",
# ...
#     "ULDR_parents_for_vMF_spherical",
#     "ULDR_children_young2_for_vMF_spherical",
#     "ULDR_children_old2_for_vMF_spherical"
# )
# 2025-07-09
# Caroline Haslebacher
# this script is a collection of visualisation code snippets for this project.
