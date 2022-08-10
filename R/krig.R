
#' Make grid for kriging
#'
#' @param xrange range of xvalues
#' @param yrange range of yvalues
#' @param len length provided to expand.grid
#'
#' @return
#' @export
#'
#' @examples
range_to_grid <- function(xrange, yrange, len){
  grd <- expand.grid(x=seq(from=xrange[1], to=xrange[2], len=len),
                     y=seq(from=yrange[1], to=yrange[2], len=len))
  coordinates(grd) <- ~x+y
  gridded(grd) <- TRUE
  return(grd)
}

#' Make grid from Spatial Points Data Frame
#'
#' @param spdf Spatial Points Dataframe to create grid for kriging
#' @param n_cell number of grid cells to use when kriging
#'
#' @note code from: https://stackoverflow.com/questions/43436466/create-grid-in-r-for-kriging-in-gstat
#' @return
#' @export
#'
#' @examples
#'
spdf_to_grid <- function(spdf, n_cell = 1000){
  # make grid from spdf
  grd <- makegrid(spdf, n = n_cell)
  colnames(grd) <- c("x", "y")
  coordinates(grd) <- ~x+y

  # Next, convert the grid to `SpatialPoints` and subset these points by the polygon.
  grd_pts <- SpatialPoints(
    coords      = grd,
    proj4string = raster::crs(spdf)
  )

  # subset all points in `grd_pts` that fall within `spdf`
  krig_grd <- grd_pts[spdf, ]

  return(krig_grd)
}


#' TEMP COORDS FUNCTION
#'
#' @param coords
#' @param spdf
#' @param crop_to_spdf
#'
#' @return
#' @export
#'
#' @examples
coord_proj <- function(coords, spdf, crop_to_spdf = FALSE){
  # make df
  coords_spdf <- data.frame(x = coords[,1], y = coords[,2])

  # make into SPFD
  coordinates(coords_spdf) <- ~x+y

  # Assign CRS
  raster::crs(coords_spdf) <- raster::crs("+proj=longlat +datum=WGS84 +no_defs")

  # Transform to CRS of spdf
  coords_spdf <- spTransform(coords_spdf, raster::crs(spdf))

  # Crop points to SPDF
  if(crop_to_spdf){coords_spdf <- crop(coords_spdf, spdf)}

  return(coords_spdf)
}
