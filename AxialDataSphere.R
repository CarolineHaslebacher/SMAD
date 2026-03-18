# Graphial representations of axial data on the sphere ----

# Rotate spherical coordinates to align with a new north pole
rotate_to_oblique <- function(lon, lat, lon_pole, lat_pole) {
  # Convert degrees to radians
  deg_to_rad <- pi / 180
  rad_to_deg <- 180 / pi
  
  lon <- lon * deg_to_rad
  lat <- lat * deg_to_rad
  lon_pole <- lon_pole * deg_to_rad
  lat_pole <- lat_pole * deg_to_rad
  
  # Compute rotation matrix components
  cos_lat_pole <- cos(lat_pole)
  sin_lat_pole <- sin(lat_pole)
  
  # Cartesian coordinates of the input points
  x <- cos(lat) * cos(lon)
  y <- cos(lat) * sin(lon)
  z <- sin(lat)
  
  # Rotate the points (this is essentially T_alpha, rotation around y-axis)
  x_rot <- cos_lat_pole * x + sin_lat_pole * z
  y_rot <- y
  z_rot <- -sin_lat_pole * x + cos_lat_pole * z
  
  # Convert back to spherical coordinates
  lon_rot <- atan2(y_rot, x_rot) * rad_to_deg
  lat_rot <- asin(z_rot) * rad_to_deg
  
  return(data.frame(lon = lon_rot, lat = lat_rot))
}

# Function for stereographic projection with oblique coordinates
equirectangular_to_oblique_stereographic_points <- function(lon, lat, lon0 = 0, lat0 = 0, 
                                                            lon_pole = 0, lat_pole = 90, radius = 1) {
  # Rotate coordinates to align the new pole with the z-axis
  rotated <- rotate_to_oblique(lon, lat, lon_pole, lat_pole)
  
  # Extract rotated coordinates
  lon_rot <- rotated$lon
  lat_rot <- rotated$lat
  
  # Convert degrees to radians
  deg_to_rad <- pi / 180
  lon_rot <- lon_rot * deg_to_rad
  lat_rot <- lat_rot * deg_to_rad
  lon0 <- lon0 * deg_to_rad
  lat0 <- lat0 * deg_to_rad
  
  # Compute common terms for stereographic projection
  cos_c <- sin(lat0) * sin(lat_rot) + cos(lat0) * cos(lat_rot) * cos(lon_rot - lon0)
  
  # Avoid division by zero (points not visible)
  denominator <- 1 + cos_c
  visible <- denominator > 0
  
  # Compute projected X and Y
  X <- 2 * radius * cos(lat_rot) * sin(lon_rot - lon0) / denominator
  Y <- 2 * radius * (cos(lat0) * sin(lat_rot) - sin(lat0) * cos(lat_rot) * cos(lon_rot - lon0)) / denominator
  
  # Mask points not visible
  X[!visible] <- NA
  Y[!visible] <- NA
  
  return(data.frame(X = X, Y = Y, visible = visible))
}

# Updated orthographic projection function for oblique projection
equirectangular_to_oblique_orthographic_points <- function(lon, lat, lon0 = 0, lat0 = 0, 
                                                           lon_pole = 0, lat_pole = 90, radius = 1) {
  # Rotate coordinates to align the new pole with the z-axis
  # NOTE: this is not used anymore, since I realised that ShowAxialDataSphere always keeps the north pole up.
  # could be deleted.
  rotated <- rotate_to_oblique(lon, lat, lon_pole, lat_pole)
  
  # Apply the standard orthographic projection to the rotated coordinates
  lon_rot <- rotated$lon
  lat_rot <- rotated$lat
  
  # Convert degrees to radians
  deg_to_rad <- pi / 180
  lon_rot <- lon_rot * deg_to_rad
  lat_rot <- lat_rot * deg_to_rad
  lon0 <- lon0 * deg_to_rad
  lat0 <- lat0 * deg_to_rad
  
  # Orthographic projection equations
  cos_c <- sin(lat0) * sin(lat_rot) + cos(lat0) * cos(lat_rot) * cos(lon_rot - lon0)
  
  # Check visibility of points (cos_c >= 0)
  visible <- cos_c >= 0
  
  # Compute projected X and Y
  X <- radius * cos(lat_rot) * sin(lon_rot - lon0)
  Y <- radius * (cos(lat0) * sin(lat_rot) - sin(lat0) * cos(lat_rot) * cos(lon_rot - lon0))
  
  # Mask points not visible from the orthographic projection
  X[!visible] <- NA
  Y[!visible] <- NA
  
  return(data.frame(X = X, Y = Y, visible = visible))
}


back_to_raster <- function(projected_visible, x_range, y_range, output_res, values_visible){
  # Step 1: Define the grid
  # x_range <- range(projected_visible$X, na.rm = TRUE)
  # y_range <- range(projected_visible$Y, na.rm = TRUE)

  # Define resolution and output raster dimensions

  ncols_out <- ceiling((x_range[2] - x_range[1]) / output_res)
  nrows_out <- ceiling((y_range[2] - y_range[1]) / output_res)

  # Create an empty raster
  output_raster <- raster(
    xmn = x_range[1], xmx = x_range[2], 
    ymn = y_range[1], ymx = y_range[2], 
    nrows = nrows_out, ncols = ncols_out
  )

  # Initialize raster values to NA
  values(output_raster) <- NA

  # Step 2: Map each projected point to the nearest raster cell
  # Compute raster cell indices for each projected point
  # The Y-axis in a raster grid increases from bottom to top, but projected Y-coordinates in a geographic coordinate system (like latitude) typically decrease from top to bottom (e.g., from the North Pole to the South Pole).
  # As a result, the top row corresponds to the maximum Y value (y_range[2]), and the row index decreases as Y decreases.
  col_indices <- floor((projected_visible$X - x_range[1]) / output_res) + 1
  row_indices <- nrows_out - floor((projected_visible$Y - y_range[1]) / output_res)

  # Ensure indices are within valid bounds
  valid <- col_indices >= 1 & col_indices <= ncols_out & 
          row_indices >= 1 & row_indices <= nrows_out

  # Filter points within the raster bounds
  col_indices <- col_indices[valid]
  row_indices <- row_indices[valid]
  values_visible <- values_visible[valid]

  # Convert row/col indices to linear indices
  linear_indices <- cellFromRowCol(output_raster, row_indices, col_indices)

  # Step 3: Assign the values to the raster
  # Handle duplicate assignments (e.g., average or overwrite)
  if (any(duplicated(linear_indices))) {
    # Aggregate values if duplicates exist (e.g., use mean)
    aggregated_values <- tapply(values_visible, linear_indices, mean)
    values(output_raster)[as.numeric(names(aggregated_values))] <- aggregated_values
  } else {
    # Assign directly if no duplicates
    values(output_raster)[linear_indices] <- values_visible
  }
  return(output_raster)
}

#%%

equirectangular_to_orthographic_raster_indices <- function(raster, lon0 = 0, lat0 = 0, lon_pole = 0, lat_pole = 0, radius = 1, output_res = 50000, projection_type = "orthographic", wsize=1) {
  # for debudding
    #   lon0 = 0 
    #   lat0 = 0
    #   radius = 2632344.9707
    #   output_res = 10000
    # projection_type <- "stereographic" # or "orthographic"  # or 
  # Get raster dimensions
  nrows <- raster::nrow(raster)
  ncols <- raster::ncol(raster)
  
  # Compute geographic coordinates from raster indices
  extent_raster <- raster::extent(raster)
  lon_step <- (extent_raster[2] - extent_raster[1]) / ncols  # Step size for longitude
  lat_step <- (extent_raster[4] - extent_raster[3]) / nrows  # Step size for latitude
  
  # Indices for each cell in the raster
  # lon and lat should be in degrees
  lon_indices <- seq(1, ncols)
  lat_indices <- seq(1, nrows)
  
  # Convert indices to geographic coordinates
  lon <- extent_raster[1] + (lon_indices - 0.5) * lon_step
  lat <- extent_raster[4] - (lat_indices - 0.5) * lat_step
  
  # Create full grid of geographic coordinates
  lonlat <- expand.grid(lat = lat, lon = lon)
  
  # Flatten raster values
  values_matrix <- raster::as.matrix(raster)
  values <- as.vector(values_matrix)
  
  if (projection_type == "orthographic") {
      # Project geographic coordinates to orthographic coordinates
      projected <- equirectangular_to_oblique_orthographic_points(lonlat$lon, lonlat$lat, lon0, lat0, lon_pole, lat_pole, radius)
      
  } else if (projection_type == "stereographic") {
      # calculate depending on viewing window
      # extent <- list(xmin = lon0-50, xmax = lon0+50, ymin = lat0-50, ymax = lat0+50) # xmin is lon, ymin is lat
      # # Filter points within the extent
      # valid_indices <- lonlat$lon >= extent$xmin & lonlat$lon <= extent$xmax &
      #                 lonlat$lat >= extent$ymin & lonlat$lat <= extent$ymax
      # lonlat <- lonlat[valid_indices, ]
      # values <- values[valid_indices]
      # Project points using stereographic projection
      projected <- equirectangular_to_oblique_stereographic_points(lonlat$lon, lonlat$lat, lon0, lat0, lon_pole, lat_pole, radius)
  } else {
      stop("Invalid projection type specified.")
  }

  # problems I see: the projected' is a dataframe with the lonlat grid.
  # It does not correspond to the order or the indices, does it?
  # why in the first place are we constructing a lonlat grid, but then only use lon and lat separately? --> because they are 'duplicated'
  # if they are duplicated, how are we extracting the 'real' values?
  # can we somehow add the values for each lonlat point into the expand.grid?
  # or do we simply need to swap lon and lat in expand.grid? ----> YES!! that was it.
  # Retain visible points and their values
  visible_indices <- which(projected$visible)
  projected_visible <- projected[visible_indices, ]
  values_visible <- values[visible_indices]

  # Replace NAs with a default value (e.g., 0)
  values_visible[is.na(values_visible)] <- 0

  # Interpolate onto the new orthographic grid
  x_range <- range(projected_visible$X, na.rm = TRUE)
  y_range <- range(projected_visible$Y, na.rm = TRUE)
  
  
 	# prepare for plotting:
	# Normalize the values for color mapping
	value_range <- range(values_visible, na.rm = TRUE)
	normalized_values <- (values_visible - value_range[1]) / diff(value_range)

	# Define a grayscale palette (256 shades)
	grayscale_palette <- colorRampPalette(c("black", "white"))
	colors <- grayscale_palette(256)[as.numeric(cut(normalized_values, breaks = 256))]

	# check by plotting:
	# png(paste("Eruopa_ortho_proj_lon0_", lon0, "_lat0_", lat0, "_lonpole_", lon_pole, "_latpole_", lat_pole ,".png"), width = 400, height = 400)
	# # Plot the projected points with grayscale colors representing the values
	# plot(projected_visible$X, projected_visible$Y, col = colors, pch = 20,
	#     xlab = "Projected X", ylab = "Projected Y",
	#     main = "Projected Points with Grayscale Values")
	# dev.off()

	# If the input basemap is small enough, I do not have to interpolate (which costs a lot of time!)

	# recreate a raster:
	# compute output_res with current grid (should have very similar ncol and nrow)
	ncols_out <- sqrt(length(projected_visible$X))
	nrows_out <- sqrt(length(projected_visible$Y))

  if (projection_type == "orthographic") {
	# note: we need a higher output resolution to decrease the cell size (because we do not interpolate, which takes time)
	output_res <- (x_range[2] - x_range[1]) / sqrt(length(projected_visible$X)) *1.8

  } else if (projection_type == "stereographic") {
      # define the range (can depend on window size)
    #   wsize <- 2
      x_range <- c(-wsize, wsize)
      y_range <- c(-wsize, wsize)
      # Project points using stereographic projection
      xres <- (x_range[2] - x_range[1]) / sqrt(length(projected_visible$X)) *4 
	  yres <- (y_range[2] - y_range[1]) / sqrt(length(projected_visible$Y)) *4 
	  # take mean of xres and yres to balance
	  output_res <- (xres+yres)/2 
	  # print(output_res)
	  } else {
      stop("Invalid projection type specified.")
	  }

  	output_raster <- back_to_raster(projected_visible, x_range, y_range, output_res, values_visible)

	# Step 4: Plot the raster to verify
	# plot(output_raster, col = grey.colors(256, start = 0, end = 1),  main = "Raster Without Interpolation")

  return(output_raster)
}

RAlpha <- function(z)
	# Auxiliary function.
	# For a vector z = z[1:2], this procedure returns its
	# standard Euclidean norm r and an angle alpha in [0,2*pi)
	# such that z == r * c(cos(alpha),sin(alpha)).
{
	r <- sqrt(sum(z^2))
	if (r <= 0){
		alpha <- 0
	}else{
		zn <- z/r
		if (zn[1] >= abs(zn[2])){
			alpha <- asin(zn[2])
			if (alpha < 0){
				alpha <- 2*pi + alpha
			}
		}
		if (zn[1] <= -abs(zn[2])){
			alpha <- pi - asin(zn[2])
		}
		if (zn[2] > abs(zn[1])){
			alpha <- acos(zn[1])
		}
		if (zn[2] < -abs(zn[1])){
			alpha <- 2*pi - acos(zn[1])
		}
	}
	return(c('r'=r,'alpha'=alpha))
}

Rotation <- function(x0)
	# For a unit vector x0 in dimension 3, this procedure
	# determines angles alpha, beta and orthogonal 3x3
	# matrices Talpha and Tbeta with determinant 1 such that:
	# - Talpha describes a rotation around the x3-axis
	#   by angle alpha,
	# - Tbeta describes a rotation around the x2-axis
	#   by angle beta,
	# - Tbeta %*% Talpha %*% x0 == c(1,0,0)
	# The procedure returns the angles alpha, beta and the
	# matrix B0 = Tbeta %*% Talpha.
{
	if (x0[3] <= -1){
		# x0 is the south pole.
		# Rotate sphere around the x2-axis such that the
		# south pole becomes the point c(1,0,0):
		alpha <- 0
		beta <- pi/2
		Talpha <- diag(3)
		Tbeta <- matrix(c(0,0,1,0,1,0,-1,0,0),3,3)
	}else{
		if (x0[3] < 1){
			# x0 is neither the south nor the north pole.
			# Rotate sphere around the x3 axis such that
			# x0[2] == 0:
			alpha <- -RAlpha(x0[1:2])[2]
			Talpha <- matrix(0,3,3)
			Talpha[3,3] <- 1
			Talpha[1,1] <- cos(alpha)
			Talpha[2,1] <- sin(alpha)
			Talpha[1,2] <- -sin(alpha)
			Talpha[2,2] <- cos(alpha)
			# Rotate sphere around the x2 axis such that
			# x0[1] == 1:
			beta <- -asin(x0[3])
			Tbeta <- matrix(0,3,3)
			Tbeta[2,2] <- 1
			Tbeta[1,1] <- cos(beta)
			Tbeta[3,1] <- sin(beta)
			Tbeta[1,3] <- -sin(beta)
			Tbeta[3,3] <- cos(beta)
			Tba <- Tbeta %*% Talpha
		}else{
			# x0 is the north pole.
			# Rotate sphere around the x2-axis such that the
			# north pole becomes the point c(1,0,0):
			alpha <- 0
			Talpha <- diag(3)
			beta <- -pi/2
			Tbeta <- matrix(c(0,0,-1,0,1,0,1,0,0),3,3)
		}
	}
	return(list(alpha=alpha,beta=beta,B0 = Tbeta %*% Talpha))
}


# Graphial representations of axial data on the sphere ----

# CH added these for plotting the stresses (custom PDFs) on the sphere.
# NOTE: these functions use spherical coordinates based on the parametrisation
# X1 = np.cos(lat)*np.cos(lon)
# X2 = np.cos(lat)*np.sin(lon)
# X3 = np.sin(lat)
# with lat e (-pi/2, pi/2]
# lon e [0, 2pi)
#
project_north <- function(theta, phi) {
  # theta : longitude
  # phi : latitude
  # Project (0, 0, 1) onto the tangential plane defined by lon (theta), lat (phi)
  
  N1 <- -sin(phi) * cos(theta)
  N2 <- -sin(phi) * sin(theta)
  N3 <- cos(phi)
  
  return(c(N1, N2, N3))
}

construct_east <- function(theta, phi) {
  # theta : longitude
  # phi : latitude
  # Construct the east direction by using the normal vector of the tangent and the north projected vector.
  
  E1 <- -sin(theta)
  E2 <- cos(theta)
  E3 <- 0 * theta # We multiply with theta in case it is an array as input
  
  return(c(E1, E2, E3))
}

alpha_to_vec <- function(alpha, lon, lat) {
  # alpha : angle defined in clockwise direction from North
  # lon : longitude
  # lat : latitude
  # output : a normalized vector of the direction that corresponds to angle alpha in the tangential plane
  
  project_north_vec <- project_north(lon, lat)
  construct_east_vec <- construct_east(lon, lat)
  
	# vec <- cos(alpha) * project_north_vec + sin(alpha) * construct_east_vec
	# for this code that uses Euclidean vectors, we simply need to change cos and sin
	vec <- sin(alpha) * project_north_vec + cos(alpha) * construct_east_vec
  absl <- sqrt(sum(vec^2)) # normalization
  
  return(vec / absl) # Normalize the vector by dividing by its length
}

# Create a function wrapper to apply element-wise
apply_alpha_to_vec <- function(a, l, b) {
  alpha_to_vec(a, l, b)
}

# changes: plotting falph (which defines the PDF) on the frontside
PlotAxialDataSphere <- function(X,V=NULL,
								scale=1,
								xysize=1.1,
								vfactor=0.1,
								alpha=0,
								beta=0,
								Tba=diag(3),
								Xcolours=c('black','gray'),
								Vcolours=c('forestgreen','gray'),
								Polecolours=c('black','gray'),
								showEquator=TRUE,
								showPoles=TRUE,
								showAxes=TRUE,
								basemapImage = NULL,
								Frame=NULL,
								Fcolours=c('blue','gray'),
								falph=NULL,
								ffactor=1,
								lwdVfactor=1,
								Ccolours=c('yellow', 'gray'))
								# note: basemapImage must be a raster()

# Auxiliary function for ShowAxialDataSphere()
{
	# Extra emphasis for V-lines in the very front:
	# lwdVfactor <- 1
	# 0: no extra emphasis, lwd=1
	# L: lwd ranges from 1 to 1+L
	if (showEquator){
		phi0 <- seq(0,2*pi,length.out=201)
		Equ0 <- cbind(cos(phi0),sin(phi0),0)
		phiNS <- seq(0,2*pi,length.out=ceiling(200*sqrt(.75)))
		EquN <- cbind(sqrt(.75)*cos(phiNS),sqrt(.75)*sin(phiNS),0.5)
		EquS <- cbind(sqrt(.75)*cos(phiNS),sqrt(.75)*sin(phiNS),-0.5)
	}
	if (showPoles){
		NP <- c(0,0,1)
	}

	Xtmp <- X %*% t(Tba)
	if (!is.null(V)){
		Vtmp <- V %*% t(Tba)
	}
	xylim <- scale*xysize*c(-1,1)
	if (showEquator){
		Equ0tmp <- Equ0 %*% t(Tba)
		plot(Equ0tmp[,2],Equ0tmp[,3],
			 pch='.',col=Xcolours[1 + (Equ0tmp[,1] < 0)],
			 xlab='',ylab='',
			 xlim=xylim,xaxs='i',
			 ylim=xylim,yaxs='i')
		EquNtmp <- EquN %*% t(Tba)
		points(EquNtmp[,2],EquNtmp[,3],
			   pch='.',col=Xcolours[1 + (EquNtmp[,1] < 0)])
		EquStmp <- EquS %*% t(Tba)
		points(EquStmp[,2],EquStmp[,3],
			   pch='.',col=Xcolours[1 + (EquStmp[,1] < 0)])
	}else{
		plot(0,0,col='white',
			 xlab='',ylab='',
			 xlim=xylim,xaxs='i',
			 ylim=yylim,yaxs='i')
	}
	# Plot data on backside:
	backside <- which(Xtmp[,1] < 0)
	points(Xtmp[backside,2],Xtmp[backside,3],pch=16,
		   col=Xcolours[2])
	if (!is.null(V)){
		for (i in backside){
			tmpx <- Xtmp[i,2] + vfactor*scale*Vtmp[i,2]*c(-1,1)
			tmpy <- Xtmp[i,3] + vfactor*scale*Vtmp[i,3]*c(-1,1)
			lines(tmpx,tmpy,col=Vcolours[2])
            # CH: I think here I can add polygons
		}
	}
	if (showPoles){
		NPtmp <- Tba %*% NP
		if (NPtmp[1] < 0){
			text(NPtmp[2],NPtmp[3],'N',
				 col=Polecolours[2])
		}
		if (NPtmp[1] > 0){
			text(-NPtmp[2],-NPtmp[3],'S',
				 col=Polecolours[2])
			
		}
	}
	if (showAxes){
		# Plot axes:
		for (j in 1:3){
			lines(c(0,-Tba[2,j]),c(0,-Tba[3,j]),
				  lty=3,col='gray')
			lines(c(0,Tba[2,j]),c(0,Tba[3,j]),
				  lty=1,col='gray')
		}
	}
	# Plot the basemap as the first layer
	if (!is.null(basemapImage)) {		
		# Note: radius = 1 works perfectly fine
		output_raster <- equirectangular_to_orthographic_raster_indices(basemapImage, lon0 = -alpha*180/pi-180, lat0 = -beta*180/pi, lon_pole = 0, lat_pole = 0, radius = 1)
		# convert output_raster to a matrix for plotting
		output_basemap <- raster::as.matrix(output_raster) 
		# normalise basemapImage (we need this...)
		#  could we normalize to full input? basemapImage_raster <- raster::as.matrix(basemapImage) 
		output_basemap <- (output_basemap - min(output_basemap, na.rm=TRUE))/(max(output_basemap, na.rm=TRUE) - min(output_basemap, na.rm=TRUE))
		rasterImage(output_basemap, -1, -1, 1, 1)
	}
	# Plot data on frontside:
	frontside <- which(Xtmp[,1] >= 0)
	points(Xtmp[frontside,2],Xtmp[frontside,3],pch=16,
		   col=Xcolours[1])
	if (!is.null(V)){
		for (i in frontside){
			if (!is.null(falph)){ 
				#  & !is.null(Valphas) Valphas are not needed any longer
				# reconstruct original lat and lon for getting Valphas (which are just the transformed alpha values so they align again with falph)
				latorig <- asin(Xtmp[i,3]) # keep it in radian, not deg. 180/pi*
				lonorig <- atan2(Xtmp[i,2], Xtmp[i,1])
				# I figured it is wrong to backtransform to Geo-coordinates.
				# calculate Valphas using lon and lat of the current display location:
				nalphas <- length(falph[i,])
				alphas <- 2*pi*(0:(nalphas-1))/nalphas # e.g. 2*pi*(0:999)/1000
				lon <- matrix(lonorig, length(alphas))
				lat <- matrix(latorig, length(alphas))
				Valphas <- t(mapply(apply_alpha_to_vec, alphas, lon, lat)) 
				# print(Valphas%*%Xtmp[i,]) # all elements need to be zero. do not take the sum, as the elements will cancel out in each case.
				# for each point, we have on Valpha and one falph
				# the Valpha is simply the Vtmp
				polygon(x=Xtmp[i,2] + vfactor*ffactor*sqrt(falph[i,])*scale*Valphas[,2],
						y=Xtmp[i,3] + vfactor*ffactor*sqrt(falph[i,])*scale*Valphas[,3], 
				col=Ccolours[1])
			}
			# we also add lines on top (for example from observations), which creates an overleay
			tmpx <- Xtmp[i,2] +
				vfactor*scale*Vtmp[i,2]*c(-1,1)
			tmpy <- Xtmp[i,3] +
				vfactor*scale*Vtmp[i,3]*c(-1,1)
			lines(tmpx,tmpy,col=Vcolours[1],
				  lwd=1+lwdVfactor*Xtmp[i,1])
		}
	}
	if (showPoles){
		if (NPtmp[1] >= 0){
			text(NPtmp[2],NPtmp[3],'N',
				 col=Polecolours[1])
		}
		if (NPtmp[1] <= 0){
			text(-NPtmp[2],-NPtmp[3],'S',
				 col=Polecolours[1])
		}
	}
	# Show Frame:
	if (!is.null(Frame)){
		lines(Frame[,2],Frame[,3],
			  col=Fcolours[1 + (Frame[,1] < 0)])
	}
}


# changes on 2025-06-09: included falph and Ccolours which are passed on to PlotAxialDataSphere
ShowAxialDataSphere <- function(X,V=NULL,x0=c(1,0,0),
								zoomfactor=2^(1/4),vfactor=0.1,
								dangle=2*pi/48,
								Xcolours=c('black','gray'),
								Vcolours=c('forestgreen','gray'),
								Polecolours=c('black','gray'),
								showEquator=TRUE,
								showPoles=TRUE,
								showAxes=TRUE,
								basemapImage = NULL,								
								falph=NULL,
								ffactor=1,
								lwdVfactor=1,
								Ccolours=c('yellow', 'gray'),
								use_locator=TRUE # by default, this is interactive and a user can define a locator in the window. FALSE is used to save plots easily (without the interactive 'zoom', 'pick x0' etc windows)
								)
# Main function to visualize and explore axial data
# on a sphere interactively. Eventually, the procedure
# returns the reference point x0, which is the point on
# the sphere closest to the observer.
{
	if (x0[1] == 1){
		alpha <- 0
		beta <- 0
	}else{
		tmp <- Rotation(x0)
		alpha <- tmp$alpha
		beta <- tmp$beta
	}
	origin <- matrix(0,1,2)
	scale <- 1
	# Define boundaries for interaction windows
	start.window <- cbind(c(-1.3,-1,-1,-1.3,-1.3),
						  c(1,1,1.3,1.3,1))
	zi.window <- cbind(c(1,1.3,1.3,1,1),
					   c(1,1,1.3,1.3,1))
	zo.window <- cbind(c(-1.3,-1,-1,-1.3,-1.3),
					   c(-1.3,-1.3,-1,-1,-1.3))
	end.window <- cbind(c(1,1.3,1.3,1,1),
						c(-1.3,-1.3,-1,-1,-1.3))
	xysize <- 1.3
	
	if (use_locator){

	answer <- TRUE
	# Note: answer is assigned 'False' only when we click on the 'End' button 
	while (answer){
		Talpha <- matrix(c(cos(alpha),sin(alpha),0,
						   -sin(alpha),cos(alpha),0,
						   0,0,1),3,3)
		Tbeta <- matrix(c(cos(beta),0,sin(beta),
						  0,1,0,
						  -sin(beta),0,cos(beta)),3,3)
		Tba <- Tbeta %*% Talpha
		PlotAxialDataSphere(X,V,scale,xysize,vfactor,
							alpha=alpha, beta=beta,
							Tba=Tba,
							Xcolours=Xcolours,Vcolours=Vcolours,
							Polecolours=Polecolours,
							showEquator=showEquator,showPoles=showPoles,
							showAxes=showAxes, basemapImage=basemapImage,
							falph=falph, ffactor=ffactor, lwdVfactor=lwdVfactor, Ccolours=Ccolours)

		# Instructions for next step:
		polygon(scale*start.window[,1],
				scale*start.window[,2],
				col='white',border='gray')
		text(scale*(-1.15),scale*1.15,'Pick x0')
		polygon(scale*zi.window[,1],
				scale*zi.window[,2],
				col='white',border='gray')
		text(scale*1.15,scale*1.15,'Z+')
		polygon(scale*zo.window[,1],
				scale*zo.window[,2],
				col='white',border='gray')
		text(scale*(-1.15),scale*(-1.15),'Z-')
		polygon(scale*end.window[,1],
				scale*end.window[,2],
				col='white',border='gray')
		text(scale*1.15,scale*(-1.15),'End')
		# read in the locator coordinates (x/y coords (because plot is in 2D))
		tmp <- locator(1)
		x <- tmp$x/scale
		y <- tmp$y/scale
		if (abs(x) > 1 & abs(y) > 1){
		  # 'Pick x0' field
			if (x < -1 & y > 1){
				# Pick a new reference point x0:
				tmp2 <- locator(1)
				x0new <- c(tmp2$x,tmp2$y)
				r2 <- sum(x0new^2)
				if (r2 >= 0.9801){
					x0new <- 0.99*x0new/sqrt(r2)
					r2 <- 0.9801
				}
				x0new <- c(sqrt(1 - r2),x0new)
				x0new <- t(Tba) %*% x0new
				if (x0new[3] <= -1){
					# x0new is the south pole.
					alpha <- 0
					beta <- pi/2
				}else{
					if (x0new[3] < 1){
						alpha <- -RAlpha(x0new[1:2])[2]
						beta <- -asin(x0new[3])
					}else{
						# x0new is the north pole
						alpha <- 0
						beta <- -pi/2
					}
				}
			}
		  # in the 'Z+' zoom field:
			if (x > 1 & y > 1){
				# Zoom in:
				scale <- scale/zoomfactor
				dangle <- dangle/zoomfactor
			}
			if (x < -1 & y < -1){
				# Zoom out:
				scale <- scale*zoomfactor
				dangle <- dangle*zoomfactor
			}
			if (x > 1 & y < -1){
				# End of journey.
				answer <- FALSE
			}
		}else{
			if (abs(x) > abs(y)){
				# Turn globe around north-south axis:
				alpha <- alpha + sign(x)*dangle
			}
			if (abs(y) > abs(x)){
				# Tilt north-south axis:
				beta <- beta + sign(y)*dangle
				beta <- max(min(pi/2,beta),-pi/2)
			}
		}
	}
	}

	xysize <- 1.1
	if (use_locator==FALSE){
		# calculate Tba
		Talpha <- matrix(c(cos(alpha),sin(alpha),0,
							-sin(alpha),cos(alpha),0,
							0,0,1),3,3)
		Tbeta <- matrix(c(cos(beta),0,sin(beta),
							0,1,0,
							-sin(beta),0,cos(beta)),3,3)
		Tba <- Tbeta %*% Talpha
		# adjust scale with zoomfactor
		scale <- scale/zoomfactor
	}
	PlotAxialDataSphere(X,V,scale,xysize,vfactor,
						alpha=alpha, beta=beta,
						Tba=Tba,
						Xcolours=Xcolours,Vcolours=Vcolours,
						Polecolours=Polecolours,
						showEquator=showEquator,showPoles=showPoles,
						showAxes=showAxes, basemapImage=basemapImage,
						falph=falph, ffactor=ffactor, lwdVfactor=lwdVfactor, Ccolours=Ccolours)
	
	x0 <- t(Tba)[,1]
	return(x0)
}


PlotWeightedAxialDataSphere <- function(X,V=NULL,
								N=100,Tba=diag(3),
								offset=0.07,
								scale=1,
								xysize=1.1,
								vfactor=0.1,
								Polecolours=c('black','gray'),
								showEquator=TRUE,
								showPoles=TRUE,
								showAxes=TRUE,
								alpha=0,
								beta=0,
								basemapImage=NULL,
								customweights=NULL)
# Auxiliary function for ShowWeightedAxialDataSphere()
{
	# Extra emphasis for V-lines in the very front:
	lwdVfactor <- 1
	# 0: no extra emphasis, lwd=1
	# L: lwd ranges from 1 to 1+L
	if (showEquator){
		phi0 <- seq(0,2*pi,length.out=201)
		Equ0 <- cbind(cos(phi0),sin(phi0),0)
		phiNS <- seq(0,2*pi,length.out=ceiling(200*sqrt(.75)))
		EquN <- cbind(sqrt(.75)*cos(phiNS),sqrt(.75)*sin(phiNS),0.5)
		EquS <- cbind(sqrt(.75)*cos(phiNS),sqrt(.75)*sin(phiNS),-0.5)
	}
	if (showPoles){
		NP <- c(0,0,1)
	}
	Xtmp <- X %*% t(Tba)
	if (!is.null(V)){
		Vtmp <- V %*% t(Tba)
	}
	xylim <- scale*xysize*c(-1,1)
	if (showEquator){
		Equ0tmp <- Equ0 %*% t(Tba)
		plot(Equ0tmp[,2],Equ0tmp[,3],
			 pch='.',col=Polecolours[1 + (Equ0tmp[,1] < 0)],
			 xlab='',ylab='',
			 xlim=xylim,xaxs='i',
			 ylim=xylim,yaxs='i')
		EquNtmp <- EquN %*% t(Tba)
		points(EquNtmp[,2],EquNtmp[,3],
			   pch='.',col=Polecolours[1 + (EquNtmp[,1] < 0)])
		EquStmp <- EquS %*% t(Tba)
		points(EquStmp[,2],EquStmp[,3],
			   pch='.',col=Polecolours[1 + (EquStmp[,1] < 0)])
	}else{
		plot(0,0,col='white',
			 xlab='',ylab='',
			 xlim=xylim,xaxs='i',
			 ylim=yylim,yaxs='i')
	}
	if (showPoles){
		NPtmp <- Tba %*% NP
		if (NPtmp[1] < 0){
			text(NPtmp[2],NPtmp[3],'N',
				 col=Polecolours[2])
			text(-NPtmp[2],-NPtmp[3],'S',
				 col=Polecolours[1])
		}else{
			text(-NPtmp[2],-NPtmp[3],'S',
				 col=Polecolours[2])
			text(NPtmp[2],NPtmp[3],'N',
				 col=Polecolours[1])
		}
	}
	if (showAxes){
		# Plot axes:
		for (j in 1:3){
			lines(c(0,-Tba[2,j]),c(0,-Tba[3,j]),
				  lty=3,col='gray')
			lines(c(0,Tba[2,j]),c(0,Tba[3,j]),
				  lty=1,col='gray')
		}
	}
	# Determine the weights:
	n <- dim(X)[1]
	W <- rep(0,n)
	tmp <- (Xtmp[,1] > -1)
	if (sum(tmp) <= N){
		W[tmp] <- 1
	}else{
		S <- rowSums((Xtmp[tmp,2:3]*2/(1 + Xtmp[,1]))^2)
		W[tmp] <- SetWeights(S,N)
	}
	if (!is.null(customweights)){
		W[tmp] <- W[tmp]*customweights[tmp]
		# normalize to (0,1) for display
		W[tmp] <- W[tmp]/max(W[tmp])
	}
	# Plot the basemap as the first layer
	if (!is.null(basemapImage)) {		
		# Note: radius = 1 works perfectly fine
		output_raster <- equirectangular_to_orthographic_raster_indices(basemapImage, lon0 = -alpha*180/pi-180, lat0 = -beta*180/pi, lon_pole = 0, lat_pole = 0, radius = 1)
		# convert output_raster to a matrix for plotting
		output_basemap <- raster::as.matrix(output_raster) 
		# normalise basemapImage (we need this...)
		#  could we normalize to full input? basemapImage_raster <- raster::as.matrix(basemapImage) 
		output_basemap <- (output_basemap - min(output_basemap, na.rm=TRUE))/(max(output_basemap, na.rm=TRUE) - min(output_basemap, na.rm=TRUE))
		rasterImage(output_basemap, -1, -1, 1, 1)
	}
	# Plot data, starting in the back:
	colXYV=gray((1-offset)*(1-W))
	ord <- order(W)
	for (i in 1:n){
		# only plot if not on backside:
		if (Xtmp[ord[i],1] > 0 ){
		points(Xtmp[ord[i],2],Xtmp[ord[i],3],pch=16,
			   col=colXYV[ord[i]])
		if (!is.null(V)){
			tmpx <- Xtmp[ord[i],2] +
				vfactor*scale*Vtmp[ord[i],2]*c(-1,1)
			tmpy <- Xtmp[ord[i],3] +
				vfactor*scale*Vtmp[ord[i],3]*c(-1,1)
			lines(tmpx,tmpy,
				  lwd=1+lwdVfactor*max(Xtmp[ord[i],1],0),
				  col=colXYV[ord[i]])
		}
		}
	}
}


ShowWeightedAxialDataSphere <- function(X,V=NULL,x0=c(1,0,0),
	N=100,
	offset=0.07,
	zoomfactor=2^(1/4),vfactor=0.1,
	dangle=2*pi/48,
	Polecolours=c('black','gray'),
	showEquator=TRUE,
	showPoles=TRUE,
	showAxes=TRUE,
	basemapImage=NULL,
	customweights=NULL,
	use_locator=TRUE) # by default, this is interactive and a user can define a locator in the window. FALSE is used to save plots easily (without the interactive 'zoom', 'pick x0' etc windows))
# Similar functionality as ShowAxialDataSphere, but
# now focussing on the local weights for the current
# reference point x0.
{
	if (x0[1] == 1){
		alpha <- 0
		beta <- 0
	}else{
		tmp <- Rotation(x0)
		alpha <- tmp$alpha
		beta <- tmp$beta
	}
	# not there in latest version: origin <- matrix(0,1,2)
	scale <- 1
	start.window <- cbind(c(-1.3,-1,-1,-1.3,-1.3),
						  c(1,1,1.3,1.3,1))
	zi.window <- cbind(c(1,1.3,1.3,1,1),
					   c(1,1,1.3,1.3,1))
	zo.window <- cbind(c(-1.3,-1,-1,-1.3,-1.3),
					   c(-1.3,-1.3,-1,-1,-1.3))
	end.window <- cbind(c(1,1.3,1.3,1,1),
						c(-1.3,-1.3,-1,-1,-1.3))
	xysize <- 1.3

	if (use_locator){

	answer <- TRUE
	while (answer){
		Talpha <- matrix(c(cos(alpha),sin(alpha),0,
						   -sin(alpha),cos(alpha),0,
						   0,0,1),3,3)
		Tbeta <- matrix(c(cos(beta),0,sin(beta),
						  0,1,0,
						  -sin(beta),0,cos(beta)),3,3)
		Tba <- Tbeta %*% Talpha
		PlotWeightedAxialDataSphere(X,V,N,Tba,offset,
							scale,xysize,vfactor,
							Polecolours,
							showEquator,showPoles,showAxes, 
							alpha, beta, basemapImage, customweights)
		
		# Instructions for next step:
		polygon(scale*start.window[,1],
				scale*start.window[,2],
				col='white',border='gray')
		text(scale*(-1.15),scale*1.15,'Pick x0')
		polygon(scale*zi.window[,1],
				scale*zi.window[,2],
				col='white',border='gray')
		text(scale*1.15,scale*1.15,'Z+')
		polygon(scale*zo.window[,1],
				scale*zo.window[,2],
				col='white',border='gray')
		text(scale*(-1.15),scale*(-1.15),'Z-')
		polygon(scale*end.window[,1],
				scale*end.window[,2],
				col='white',border='gray')
		text(scale*1.15,scale*(-1.15),'End')
		tmp <- locator(1)
		x <- tmp$x/scale
		y <- tmp$y/scale
		if (abs(x) > 1 & abs(y) > 1){
			if (x < -1 & y > 1){
				tmp2 <- locator(1)
				x0new <- c(tmp2$x,tmp2$y)
				r2 <- sum(x0new^2)
				if (r2 >= 0.9801){
					x0new <- 0.99*x0new/sqrt(r2)
					r2 <- 0.9801
				}
				x0new <- c(sqrt(1 - r2),x0new)
				x0new <- t(Tba) %*% x0new
				if (abs(x0new[3] + 1) < 10^(-5)){
					# x0new is close to the south pole.
					alpha <- 0
					beta <- pi/2
				}else{
					if (abs(x0new[3] - 1) >= 10^(-5)){
						alpha <- -RAlpha(x0new[1:2])[2]
						beta <- -asin(x0new[3])
					}else{
						beta <- -pi/2
					}
				}
			}
			if (x > 1 & y > 1){
				scale <- scale/zoomfactor
				dangle <- dangle/zoomfactor
			}
			if (x < -1 & y < -1){
				scale <- scale*zoomfactor
				dangle <- dangle*zoomfactor
			}
			if (x > 1 & y < -1){
				answer <- FALSE
			}
		}else{
			if (abs(x) > abs(y)){
				alpha <- alpha + sign(x)*dangle
			}
			if (abs(y) > abs(x)){
				beta <- beta + sign(y)*dangle
				beta <- max(min(pi/2,beta),-pi/2)
			}
		}
	}
	}
	
	xysize <- 1.1
	if (use_locator==FALSE){
		# calculate Tba
		Talpha <- matrix(c(cos(alpha),sin(alpha),0,
						   -sin(alpha),cos(alpha),0,
						   0,0,1),3,3)
		Tbeta <- matrix(c(cos(beta),0,sin(beta),
						  0,1,0,
						  -sin(beta),0,cos(beta)),3,3)
		Tba <- Tbeta %*% Talpha
		# and adjust zoomfactor/scale
		scale <- scale/zoomfactor
	}
	PlotWeightedAxialDataSphere(X,V,N,Tba,offset,
								scale,xysize,vfactor,
								Polecolours,
								showEquator,showPoles,showAxes,
								alpha, beta, basemapImage, customweights)
	x0 <- t(Tba)[,1]
	return(x0)
}


# Stereographic projections of axial data on the sphere ----

StereogrProj <- function(x0,X,V=NULL,
						 plotdata=TRUE,
						 xysize=2,
						 Xcolour='black',
						 vfactor=0.1,
						 Vcolour='forestgreen',
						 basemapImage = NULL)
						# note: basemapImage must be a raster()
	# For axial data on the sphere (given by X,V),
	# this procedure provides their stereopgraphic projection
	# onto the tangent plane of the point x0.
	# If x0 is not the north or south pole, the sphere is
	# rotated along the pole axis such that x0 lies
	# in the x1-x3-plane with positive first component.
	# Then the sphere is rotated around the x2-axis such that
	# x0 becomes the point c(1,0,0).
	# The rows of V need not be unit vectors. They are
	# mapped onto vectors of the same length, provided
	# that V[i,] is always perpendicular to X[i,].
{
	# Determine rotation matrix such that x0 becomes c(1,0,0):
	B0 <- Rotation(x0)$B0
	Xnew <- X %*% t(B0)
	tmp <- (Xnew[,1] > -1)
	if (sum(tmp) < 2){
		print('Oops, most points are antipodes of x0!')
		return(NULL)
	}
	Xnew <- Xnew[tmp,]
	if (!is.null(V)){
		Vnew <- V[tmp,] %*% t(B0)
	}else{
		Vnew <- NULL
	}

	# calculate lon/lat (alpha and beta) for basemap
	if (x0[3] <= -1){
		alpha <- 0
		beta <- pi/2
	}else{
		if (x0[3] < 1){
			alpha <- -RAlpha(x0[1:2])[2]
			beta <- -asin(x0[3])
		}else{
			alpha <- 0
			beta <- -pi/2
		}
	}

	# Now comes the stereographic projection onto the
	# tangent plane of c(1,0,0):
	lX <- 2/(1 + Xnew[,1])
	X2 <- lX*Xnew[,2:3]
	if (!is.null(V)){
		V2 <- Vnew[,2:3] -
			0.5*lX*Vnew[,1]*Xnew[,2:3]
	}else{
		V2 <- NULL
	}

	# Plot of the projected data:
	if (plotdata){
		plot(xysize*c(-1,1),xysize*c(-1,1),col='white',
			 xlab='',ylab='',
			 xaxs='i',yaxs='i')
			# Plot the basemap as the first layer
			if (!is.null(basemapImage)) {		
				# print(alpha)
				# print(beta)
				# Note: radius = 1 works perfectly fine
				# same function as orthographic, but with projection_type = "stereographic"
				output_raster <- equirectangular_to_orthographic_raster_indices(basemapImage, lon0 = -alpha*180/pi-180, lat0 = -beta*180/pi, lon_pole = 0, lat_pole = 0, radius = 1, projection_type = "stereographic", wsize=xysize)
				# convert output_raster to a matrix for plotting
				# output_basemap <- raster::as.matrix(output_raster) 
				# normalise basemapImage (we need this...)
				#  could we normalize to full input? basemapImage_raster <- raster::as.matrix(basemapImage) 
				# output_basemap <- (output_basemap - min(output_basemap, na.rm=TRUE))/(max(output_basemap, na.rm=TRUE) - min(output_basemap, na.rm=TRUE))
				# rasterImage(output_basemap, -1, -1, 1, 1)
				plot(output_raster, col = grey.colors(256, start = 0, end = 1), legend = FALSE)
			}
		points(X2[,1],X2[,2],pch=16,col=Xcolour)
		if (!is.null(V2)){
			for (i in which(tmp)){
				tmpx <- X2[i,1] + vfactor*xysize*V2[i,1]*c(-1,1)
				tmpy <- X2[i,2] + vfactor*xysize*V2[i,2]*c(-1,1)
				lines(tmpx,tmpy,col=Vcolour)
			}
		}
	}
	
	return(list(X2=X2,V2=V2,B0=B0))
}

InvStereogrProj <- function(x0,X2,V2=NULL)
	# This procedure provides the inverse of the
	# stereographic projection with respect to a
	# target point x0 on the sphere.
{
	# Inverse stereographic projection of the data:
	if (is.vector(X2)){
		X2 <- matrix(X2,1,2)
	}
	mX <- 1/(1 + rowSums(X2^2)/4)
	X <- cbind(2*mX - 1,mX*X2)
	if (!is.null(V2)){
		if (is.vector(V2)){
			V2 <- matrix(V2,1,2)
		}
		XV <- rowSums(X2*V2)
		V <- cbind(-mX*XV,V2 - 0.5*mX*XV*X2)
	}else{
		V <- NULL
	}
	
	# Find a rotation that maps x0 onto c(1,0,0):
	B0 <- Rotation(x0)$B0
	X <- X %*% B0
	if (!is.null(V)){
		V <- V %*% B0
	}
	return(list(X=X,V=V))
}

StereogrTour <- function(x0=c(1,0,0),X,V=NULL,
						 zoomfactor=2^(1/4),
						 xysize=1,
						 vfactor=0.1,
						 Xcolour='black',
						 Vcolour='forestgreen',
						 ShowPoles=TRUE,
						 basemapImage=NULL)
	# Visualize and explore axial data on a sphere
	# interactively while looking at the stereographic
	# projection of the data.
{
	zi.window <- cbind(c(0.8,1,1,0.8,0.8),
					   c(0.8,0.8,1,1,0.8))
	zo.window <- cbind(c(-1,-0.8,-0.8,-1,-1),
					   c(-1,-1,-0.8,-0.8,-1))
	end.window <- cbind(c(0.8,1,1,0.8,0.8),
						c(-1,-1,-0.8,-0.8,-1))
	if (ShowPoles){
		Poles <- rbind(c(0,0,1),c(0,0,-1))
	}
	xysize <- 1
	
	answer <- TRUE
	while (answer){
		stp <- StereogrProj(x0,X,V,xysize=xysize,
					 Xcolour=Xcolour,
					 Vcolour=Vcolour,vfactor=vfactor, 
					 basemapImage=basemapImage)
		if (ShowPoles && abs(x0[3]) < 1){
			Poles2 <- StereogrProj(x0,X=Poles,plotdata=FALSE)$X2
			text(Poles2[1,1],Poles2[1,2],'N')
			text(Poles2[2,1],Poles2[2,2],'S')
		}
		if (ShowPoles && x0[3] >= 1){
			text(0,0,'N')
		}
		if (ShowPoles && x0[3] <= -1){
			text(0,0,'S')
		}
		title(main=paste('x0 = [',round(x0[1],3),
						 ',',round(x0[2],3),
						 ',',round(x0[3],3),"]'",sep=''))
		# Instructions for next step:
		polygon(xysize*zi.window[,1],
				xysize*zi.window[,2],
				col='white',border='gray')
		text(xysize*0.9,xysize*0.9,'Z+')
		polygon(xysize*zo.window[,1],
				xysize*zo.window[,2],
				col='white',border='gray')
		text(xysize*(-0.9),xysize*(-0.9),'Z-')
		polygon(xysize*end.window[,1],
				xysize*end.window[,2],
				col='white',border='gray')
		text(xysize*0.9,xysize*(-0.9),'End')
		tmp <- locator(1)
		x <- tmp$x/xysize
		y <- tmp$y/xysize
		if (abs(x) > 0.8 & abs(y) > 0.8){
			if (x < -0.8 & y > 0.8){
				# Pick a new reference point x0:
				x02 <- c(tmp$x,tmp$y)
				x0 <- InvStereogrProj(x0,X2=x02)$X
			}
			if (x > 0.8 & y > 0.8){
				# Zoom in:
				xysize <- xysize/zoomfactor
			}
			if (x < -0.9 & y < -0.8){
				# Zoom out:
				xysize <- xysize*zoomfactor
			}
			if (x > 0.8 & y < -0.8){
				# End of journey.
				answer <- FALSE
			}
		}else{
			# Pick a new reference point x0:
			x02 <- c(tmp$x,tmp$y)
			x0 <- InvStereogrProj(x0,X2=x02)$X
		}
	}
	stp <- StereogrProj(x0,X,V,xysize=xysize,
				 Xcolour=Xcolour,
				 Vcolour=Vcolour,vfactor=vfactor, 
				 basemapImage=basemapImage)
	if (ShowPoles && abs(x0[3]) < 1){
		Poles2 <- StereogrProj(x0,X=Poles,plotdata=FALSE)$X2
		text(Poles2[1,1],Poles2[1,2],'N')
		text(Poles2[2,1],Poles2[2,2],'S')
	}
	if (ShowPoles && x0[3] >= 1){
		text(0,0,'N')
	}
	if (ShowPoles && x0[3] <= -1){
		text(0,0,'S')
	}
	title(main=paste('x0 = [',round(x0[1],3),
					 ',',round(x0[2],3),
					 ',',round(x0[3],3),"]'",sep=''))
	return(x0)
}




# Smoothing of axial data on a sphere ----

GreedySubsetSphere <- function(X,size,V=NULL, returnJ=FALSE)
	# Auxiliary function.
	# For a data matrix X (n x 3) with n > 1, this procedure
	# determines stepwise a submatrix X0 = X[J,] with size rows
	# such that the single rows are approximately as far apart
	# as possible.
	# to debug:
	# X <- X0_Europa
	# size <- 150
	# returnJ <- TRUE
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
	J[1] <- sample(n,1) # random sample from n. take 1
	dist <- 2*(1 - colSums(t(X)*X[J[1],]))
	for (i in 2:size){
		J[i] <- which.max(dist)
		dist <- pmin(dist,2*(1 - colSums(t(X) * X[J[i],])))
	}
	J <- sort(J)
	if (returnJ==FALSE){
		if (is.null(V)){
			return(X[J,])
		}else{
			return(list(X0=X[J,],V0=V[J,]))
		}
	}else{
		if (is.null(V)){
			return(list(X0=X[J,], J=J))
		}else{
			return(list(X0=X[J,],V0=V[J,], J=J))
		}		
	}
}

SmoothAxialDataSphere <- function(X0,X,V,N=NULL,margin=10^(-5),
								  showsteps=TRUE, regfac=0, customweights=NULL)
{
	if (is.vector(X0)){
		X0 <- matrix(X0,1,3)
	}
	# W0const <- array(0,dim(X0))
	# M0const <- array(0,dim(X0))
	# W0lin   <- array(0,dim(X0))
	# M0lin   <- array(0,dim(X0))
	# W0quadr <- array(0,dim(X0))
	# M0quadr <- array(0,dim(X0))
	nn <- dim(X0)[1]
	W02const <- matrix(0,nn,2)
	W0const  <- matrix(0,nn,3)
	M0const  <- matrix(0,nn,3)
	W02lin   <- matrix(0,nn,2)
	W0lin    <- matrix(0,nn,3)
	M0lin    <- matrix(0,nn,3)
	W02quadr <- matrix(0,nn,2)
	W0quadr  <- matrix(0,nn,3)
	M0quadr  <- matrix(0,nn,3)
	for (j in 1:nn){
		# Stereographic projection onto tangent plane at X0[j,]:
		ii <- (X %*% X0[j,] + 1 > margin)
		# print(ii)
		res <- StereogrProj(X0[j,],X[ii,],V[ii,],plotdata=FALSE)
		x02 <- c(0,0)
		X2 <- res$X2
		V2 <- res$V2
		B0 <- res$B0
		Y2 <- cbind(V2[,1]^2 - V2[,2]^2, 2*V2[,1]*V2[,2])
		# Local constant, linear and quadratic fits:
		if(!is.null(customweights)){
			PF <- LocalPolynGLM_vMF(Y=Y2,X=X2,xx=x02,N=N,xsf=1, regfac=regfac, customweights=customweights[ii])
		}else{
			# we do not pass any customweights
			PF <- LocalPolynGLM_vMF(Y=Y2,X=X2,xx=x02,N=N,xsf=1, regfac=regfac)

		}
		tmp <- RAlpha(PF$Zcxx)
		W02const[j,] <- tmp[1]*c(cos(tmp[2]/2),sin(tmp[2]/2))
		W0const[j,] <- t(B0) %*% c(0,W02const[j,])
		tmp <- RAlpha(PF$Mcxx)
		M0const[j,] <- t(B0) %*% c(0,
			tmp[1]*c(cos(tmp[2]/2),sin(tmp[2]/2)))
		tmp <- RAlpha(PF$Zlxx)
		W02lin[j,] <- tmp[1]*c(cos(tmp[2]/2),sin(tmp[2]/2))
		W0lin[j,] <- t(B0) %*% c(0,W02lin[j,])
		tmp <- RAlpha(PF$Mlxx)
		M0lin[j,] <- t(B0) %*% c(0,
			tmp[1]*c(cos(tmp[2]/2),sin(tmp[2]/2)))
		tmp <- RAlpha(PF$Zqxx)
		W02quadr[j,] <- tmp[1]*c(cos(tmp[2]/2),sin(tmp[2]/2))
		W0quadr[j,] <- t(B0) %*% c(0,W02quadr[j,])
		tmp <- RAlpha(PF$Mqxx)
		M0quadr[j,] <- t(B0) %*% c(0,
			tmp[1]*c(cos(tmp[2]/2),sin(tmp[2]/2)))
		if (showsteps){
			print(round(c('j'=j,X0[j,],
				'Mc'=M0const[j,],
				'Ml'=M0lin[j,],
				'Mq'=M0quadr[j,]),3))
		}
	}
	return(list(X0=X0,
				W02const=W02const,
				W0const=W0const,M0const=M0const,
				W02lin=W02lin,
				W0lin=W0lin,M0lin=M0lin,
				W02quadr=W02quadr,
				W0quadr=W0quadr,M0quadr=M0quadr))
}



# Bingham distributions etc. ----

# by Caroline
PlotEmpiricalAxialAngularHistogramm <- function(ftheta,theta,
								xymax=NULL,gridsize=48)
# this function generates an angular histogram out of the given funciton ftheta with accompanying values theta 
{
	xfth <- sqrt(ftheta)*cos(theta)
	yfth <- sqrt(ftheta)*sin(theta)
	if (is.null(xymax)){
		xymax <- sqrt(max(ftheta))
	}
	xylim <- xymax*c(-1,1)
	par(cex=1.2,mai=c(0.4,0.4,0.4,0.1),mgp=c(2,0.5,0))
	plot(cos(theta),sin(theta),type='l',col='cornflowerblue',
		 xlim=xylim,ylim=xylim,
		 xaxs='i',yaxs='i',
		 xlab='',ylab='',
		 main=paste("Empirical Angular Histogram",sep=''))
	polygon(xfth,yfth,col=rgb(0.9,0.9,0.9))
	for (i in seq_along(theta)){
		lines(c(0,xfth[i]),c(0,yfth[i]),
			  col=rgb(0.5,0.5,0.5))
	}
	lines(cos(theta),sin(theta),col='cornflowerblue',lty=3)
	return()
}


AxialAngularDensity <- function(w,theta=2*pi*(0:360)/360,
								displayDensity=TRUE,
								fmax=NULL,
								gridsize=48,sun=FALSE)
{
	tmp <- RAlpha(w)
	kappa <- tmp[1]
	beta <- tmp[2]
	if(sun){
		# we set the minimum to 1 to produce weights useful to correct the sun azimuth bias
		ftheta <- exp(kappa*(1+cos(2*(theta - beta))) - lG0(kappa,d=2))
	}else{
		ftheta <- exp(kappa*cos(2*(theta - beta)) - lG0(kappa,d=2))
	}
	if (displayDensity){
		if (is.null(fmax)){
			fmax=1.05*max(ftheta)
		}
		par(cex=1.2,mai=c(0.7,0.75,0.4,0.1),mgp=c(1.7,0.5,0))
		plot(c(0,2*pi),rep(1,2),
			 type='l',col='cornflowerblue',
			 xlim=c(0,2*pi),xaxs='i',
			 ylim=c(0,fmax),yaxs='i',
			 xlab=expression(italic(theta)),
			 ylab=expression(italic(tilde(f)(theta))),
			 main=paste('w = [',
			 		   w[1],',',w[2],"]'",sep=''))
		polygon(c(theta-2*pi,theta,theta+2*pi,4*pi,-2*pi,-2*pi),
				c(ftheta,ftheta,ftheta,0,0,ftheta[1]),
				col='lightgray')
		abline(v=c(beta,beta+pi))
		lines(c(0,2*pi),rep(1,2),col='cornflowerblue',lty=3)
	}
	return(ftheta)
}

AxialAngularHistogramm <- function(w,theta=2*pi*(0:360)/360,
							xymax=NULL,gridsize=48,sun=FALSE)
{
	tmp <- RAlpha(w)
	kappa <- tmp[1]
	beta <- tmp[2]
	if(sun){
		# we set the minimum to 1 to produce weights useful to correct the sun azimuth bias
		# ftheta[ftheta<1] <- 1
		ftheta <- exp(kappa*(1+cos(2*(theta - beta))) - lG0(kappa,d=2))
	}else{
		ftheta <- exp(kappa*cos(2*(theta - beta)) - lG0(kappa,d=2))
	}
	xfth <- sqrt(ftheta)*cos(theta)
	yfth <- sqrt(ftheta)*sin(theta)
	if (is.null(xymax)){
		xymax <- sqrt(max(ftheta))
	}
	xylim <- xymax*c(-1,1)
	par(cex=1.2,mai=c(0.4,0.4,0.4,0.1),mgp=c(2,0.5,0))
	plot(cos(theta),sin(theta),type='l',col='cornflowerblue',
		 xlim=xylim,ylim=xylim,
		 xaxs='i',yaxs='i',
		 xlab='',ylab='',
		 main=paste("w = [",
		 		   w[1],',',w[2],"]'",sep=''))
	polygon(xfth,yfth,col=rgb(0.9,0.9,0.9))
	grid <- 2*pi*(1:gridsize)/gridsize
	fgrid <- exp(kappa*cos(2*(grid - beta)) - lG0(kappa,d=2))
	xfgrid <- sqrt(fgrid)*cos(grid)
	yfgrid <- sqrt(fgrid)*sin(grid)
	for (i in 1:gridsize){
		lines(c(0,xfgrid[i]),c(0,yfgrid[i]),
			  col=rgb(0.5,0.5,0.5))
	}
	lines(cos(theta),sin(theta),col='cornflowerblue',lty=3)
	return(ftheta)
}


WassersteinDiscrete <- function(p,q,CM)
	# Illustrates the optimal coupling of two
	# discrete distributions on {1,2,...,m}, given by
	# their point mass vectors p and q.
	# This procedure is not used later...
{
	m <- length(p)
	n <- length(q)
	PQ <- matrix(0,m,n)
	I <- 1
	J <- 1
	while (I <= m && J <= n){
		PQ[I,J] <- min(p[I],q[J])
		p[I] <- p[I] - PQ[I,J]
		q[J] <- q[J] - PQ[I,J]
		if (p[I] <= 0){
			I <- I+1
		}
		if (q[J] <= 0){
			J <- J+1
		}
	}
	return(PQ)
}

WassersteinAxial0 <- function(p,q,x,y)
	# Auxiliary function for WassersteinAxial()
{
	m <- length(p)
	tpi <- 2*pi
	W <- 0
	I <- 1
	J <- 1
	while (I <= m && J <= m){
		PQ <- min(p[I],q[J])
		p[I] <- p[I] - PQ
		q[J] <- q[J] - PQ
		dist <- abs(x[I] - y[J])
		W <- W + PQ*min(dist,tpi-dist)
		if (p[I] <= 0){
			I <- I+1
		}
		if (q[J] <= 0){
			J <- J+1
		}
	}
	return(W)
}

WassersteinAxial <- function(p,q)
	# For two discrete distributions P and Q on the unit circle,
	# given by point mass vectors p and q, where
	#    p[i] = P({2*pi*i/m + alpha})
	# and
	#    q[i] = Q({2*pi*i/m} + alpha})
	# for i = 1,2,...,m and an arbitrary alpha, this procedure
	# computes the Wasserstein-1-distance between P and Q with
	# respect to the usual distance along the circle.
{
	m <- length(p)
	x <- 2*pi*(1:m)/m
	y <- x
	W <- Inf
	for (k in 1:m){
		q <- c(q[-1],q[1])
		y <- c(y[-1],y[1])
		W <- min(W,WassersteinAxial0(p,q,x,y))
	}
	return(W)
}

WassersteinBingham <- function(w1,w2,m=1000,display=TRUE)
	# Computes/approximates the Wasserstein-1-distance
	# between the Bingham distributions with parameters
	# (vectors) w1 and w1.
{
	if (display){
		par(mfrow=c(2,2))
		theta <- 2*pi*(0:m)/m
		f1 <- AxialAngularDensity(w1,theta,display=TRUE)
		f2 <- AxialAngularDensity(w2,theta,display=TRUE)
		plot(theta,f1,type='l',lwd=2,
			 xlab=expression(italic(theta)),
			 ylab=expression(italic(f(theta))))
		abline(h=1,lty=3,col='blue')
		plot(theta,f2,type='l',lwd=2,
			 xlab=expression(italic(theta)),
			 ylab=expression(italic(f(theta))))
		abline(h=1,lty=3,col='blue')
		par(mfrow=c(1,1))
	}
	theta <- 2*pi*(1:m)/m
	p1 <- AxialAngularDensity(w1,theta,display=FALSE)/m
	p2 <- AxialAngularDensity(w2,theta,display=FALSE)/m
	W <- WassersteinAxial(p1,p2)
	return(W)
}

WassersteinBinghamEmpirical <- function(w,q,mult=5,display=TRUE,sun=FALSE)
	# Computes/approximates the Wasserstein-1-distance
	# between the Bingham distribution with parameter
	# (vector) w and a distribution Q on the
	# unit circle given by weights
	#    q[i] = Q((2*pi*(i - 1)/m,2*pi*i/m]), 1 <= i <= m .
{
	m <- length(q)
	q <- rep(q,each=mult)
	q <- q/sum(q) # normalization to sum(q) = 1
	mm <- m*mult
	theta <- 2*pi*((1:mm) - 0.5)/mm
	p <- AxialAngularDensity(w,theta,display=FALSE,sun=sun)
	# print(sum(p))
	# if(sun){
	# 	p[p<1] <- 1
	# 	print(sum(p))
	# 	print(max(p))
	# }
	p <- p/sum(p)
	# print(sum(p))
	W <- WassersteinAxial(p,q)
	if (display){
		ftheta <- mm*p/2/pi
		gtheta <- mm*q/2/pi
		plot(theta,gtheta,type='l',lwd=2,
			 xlab=expression(italic(theta)),
			 ylab=expression(italic(f(theta))),
			 ylim=c(0,max(ftheta,gtheta)),
			 main=paste('W1 =',round(W,4)))
		lines(theta,ftheta,col='blue')
	}
	return(W)
}
