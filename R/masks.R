
#' Create raster mask based on coordinates
#'
#' @description Creates a raster that can be used to mask areas falling outside the observation range of a dataset, as defined by coordinates and corresponding raster values
#' @param method method to create mask (can be "range", "sd", "buffer", defaults to "range"). See details for more information.
#' @return raster mask where values of 1 indicate areas that fall outside of observation range
#' @details method can either be:
#' 1. range - uses \code{\link{range_mask}}, mask all areas with values outside of the range of any of the values of the coords
#' 2. sd - uses \code{\link{sd_mask}}, mask all areas outside the mean +/- stdev*nsd of any of the values of the coords (\code{nsd} defaults to 2)
#' 3. buffer - uses \code{\link{buffer_mask}}, mask all areas outside of the buffer_width around the coords (\code{buffer_width} defaults to 0.8)
#' 4. chull - uses \code{\link{chull_mask}}, mask all areas outside a convex hull of the points
#'
#' @export
#'
#' @examples
extrap_mask <- function(coords, envlayers, method = "range", nsd = 2, buffer_width = NULL){

  if(method == "range"){
    map_mask <- range_mask(coords, envlayers)
  }

  if(method == "sd"){
    map_mask <- sd_mask(coords, envlayers, nsd = nsd)
  }

  if(method == "buffer"){
    # If no buffer width is provided, use default value of 0.8
    if(is.null(buffer_width)){buffer_width <- 0.8}
    map_mask <- buffer_mask(coords, envlayers, buffer_width = buffer_width)
  }

  if(method == "chull"){
    map_mask <- chull_mask(coords, envlayers, buffer_width = buffer_width)
  }

  # Returns raster where values of 1 are areas that should be masked, and NAs are areas that should be retained
  return(map_mask)

}

#' Create raster mask based on range of data
#' @describeIn extrap_mask mask based on range of data
#'
#' @param coords data frame of coordinates (first column should be x and second should be y)
#' @param envlayers stack of rasters with environmental values to base mask on
#'
#' @export
#'
#' @examples
range_mask <- function(coords, envlayers){
  # Extract values at all coords
  vals <- raster::extract(envlayers, coords)

  # Create empty layer of zeroes
  envmask <- envlayers*0

  # Calculate min and max values
  vals <- as.matrix(vals)
  val_max <- apply(vals, 2, max, na.rm = TRUE)
  val_min <- apply(vals, 2, min, na.rm = TRUE)

  # Loop to assign values of 1 to areas that should be masked based on the min/max vals for each layer
  for(n in 1:nlayers(envlayers)){
    envmask[[n]][envlayers[[n]] > val_max[n]] <- 1
    envmask[[n]][envlayers[[n]] < val_min[n]] <- 1
  }

  # Sum layers together to get all areas that should be masked
  map_mask <- raster::stackApply(envmask, 1, sum, na.rm=TRUE)
  # Assign values of 1 to any areas that should be masked (e.g. anything that is not 0)
  map_mask[map_mask != 0] <- 1
  # Assign NA values to any areas that should not be masked (i.e. any zeros)
  map_mask[map_mask == 0] <- NA

  return(map_mask)
}


#' Create raster mask based on mean and standard deviation of data
#' @describeIn extrap_mask mask based on mean and standard deviation of data
#'
#' @param coords data frame of coordinates (first column should be x and second should be y)
#' @param envlayers stack of rasters with environmental values to base mask on
#' @param nsd number of standard deviations to use if using the "sd" method
#'
#' @return
#' @export
#'
#' @examples
sd_mask <- function(coords, envlayers, nsd){

  vals <- raster::extract(envlayers, coords)
  # TODO [EAC]: check that this is ok
  vals <- as.matrix(vals)
  envmask <- envlayers*0

  val_mean <- apply(vals, 2, mean, na.rm = TRUE)
  val_sd <- apply(vals, 2, sd, na.rm = TRUE)

  val_max <- val_mean + val_sd*nsd
  val_min <- val_mean - val_sd*nsd

  # Loop to assign values of 1 to areas that should be masked based on the min/max vals for each layer
  for(n in 1:nlayers(envlayers)){
    envmask[[n]][envlayers[[n]] > val_max[n]] <- 1
    envmask[[n]][envlayers[[n]] < val_min[n]] <- 1
  }

  # Sum layers together to get all areas that should be masked
  map_mask <- stackApply(envmask, 1, sum, na.rm=TRUE)
  # Assign values of 1 to any areas that should be masked (i.e., anything that is not 0)
  map_mask[map_mask != 0] <- 1
  # Assign NA values to any areas that should not be masked (i.e., any zeros)
  map_mask[map_mask == 0] <- NA

  return(map_mask)

}


#' Mask rasters based on buffers around points
#' @describeIn extrap_mask mask based on buffers around points
#'
#' @param coords data frame of coordinates (first column should be x and second should be y)
#' @param envlayers stack of rasters with environmental values to base mask on
#' @param buffer_width buffer width to supply to \code{gBuffer} if using "buffer" or "chull" method. If "buffer" method is used, defaults to 0.8 and if "chull" method is used, defaults to null (no buffer)
#'
#' @return
#' @export
#'
#' @examples
buffer_mask <- function(coords, envlayers, buffer_width = 0.8){

  # Create proper coords and add projection
  colnames(coords) <- c("x", "y")
  sp::coordinates(coords) <- ~x+y
  raster::crs(coords) <- raster::crs(envlayers)

  # Add a buffer
  buff <- rgeos::gBuffer(coords, width = buffer_width)

  # Create a mask (just need one of the envlayers to do this since the values don't matter)
  map_mask <- raster::mask(envlayers[[1]], buff, inverse = TRUE)

  # Modify mask to get values of 1 for any areas that should be masked
  # (multiplying by 0 and adding 1 will make everything but the NA values equal to 1)
  map_mask <- 0*map_mask + 1

  return(map_mask)
}

#' Mask rasters based on convex hull around points
#' @describeIn extrap_mask mask based on range of data
#'
#' @param coords data frame of coordinates (first column should be x and second should be y)
#' @param envlayers stack of rasters with environmental values to base mask on
#' @param buffer_width buffer width to supply to \code{gBuffer} if using "buffer" method
#'
#' @return
#' @export
#'
#' @examples
chull_mask <- function(coords, envlayers, buffer_width = NULL){

  # Use one layer as a template
  env <- envlayers[[1]]

  # Create proper coords and add projection
  colnames(coords) <- c("x", "y")
  sp::coordinates(coords) <- ~x+y
  raster::crs(coords) <- raster::crs(env)

  # Add a buffer to coords
  if(!is.null(buffer_width)){coords <- rgeos::gBuffer(coords, width = buffer_width)}

  # Make convex hull
  chull <- rgeos::gConvexHull(coords)

  # Create mask from areas outside of hull
  env_chull_mask <- raster::mask(env, chull, inverse=TRUE)

  # Convert to ones
  map_mask <- env_chull_mask*0 + 1

  return(map_mask)

}

#' Plot mask on top of map
#' @description Plots a raster with another raster "mask" on top of it
#'
#' @param map_r raster you want masked
#' @param map_mask raster layer with 1s where you want to mask and NA everywhere else (i.e., what you want to keep, as produced by \code{\link{extrap_mask}})
#' @param RGB_cols whether the plot should be RGB-based or not
#' @param mask_col color and transparency of mask (defaults to black and alpha=0.9)
#'
#' @return plot \code{map} with areas masked based on \code{map_mask}
#' @export
#'
#' @examples
#' @seealso \code{\link{extrap_mask}}
plot_extrap_mask <- function(map_r, map_mask, RGB_cols = TRUE, mask_col = rgb(0, 0, 0, alpha = 0.9)){

  if(RGB_cols){raster::plotRGB(map_r, r = 1, g = 2, b = 3)} else {raster::plot(map_r)}

  # Plots mask as black semi-transparent layer over map
  raster::plot(map_mask, col = mask_col, add = TRUE, legend = FALSE)

}
