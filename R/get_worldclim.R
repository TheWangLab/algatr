
#' Download and merge WorldClim data for study area
#'
#' @param coords Dataframe with x and y sample coordinates.
#' @param res The resolution of WorldClim data to download; options are 0.5, 2.5, 5, and 10 arc-minutes (default = 0.5).
#' @param buff A buffer area around sample points for cropping the data layers, expressed as a proportion of the spatial extent for the coordinates (default = 0.01).
#' @param folder A directory in which to place downloaded worldclim data; if NULL then files will be downloaded to a tmp folder in the working directory (default = NULL).
#' @details
#' If res = 0.5 then the individual WorldClim tiles that cover the sample coordinates are downloaded and merged. If res > 2.5 then global layers are downloaded.
#' The buffer area maintains a large extent for the final cropped data layers around the sample coordinates. e.g. buff = 0.01 creates a 1% buffer area around the coordinates.
#' @return A SpatRaster of WorldClim layers.
#' @export
#'
#' @examples
get_worldclim <- function(coords, res = 0.5, buff = 0.01, folder = NULL){
  # Raster of worldclim tiles
  r <- raster::raster(vals = 1:60, nrows = 5, ncols = 12, ext = raster::extent(c(-180, 180, -90, 90)))

  # Make SpatialPolygons object with convex hull of coords
  ch_pts <- grDevices::chull(coords)
  ch_poly <- sp::Polygon(coords[ch_pts,])
  ch_polys <- sp::Polygons(list(ch_poly), ID = "chull")
  ch_spolys <- sp::SpatialPolygons(list(ch_polys))

  # Identify WorldClim tiles to download
  r_nums <- raster::extract(r, ch_spolys)
  r_xy <- raster::xyFromCell(r, r_nums[[1]])

  # Download and merge tiles
  message("Downloading WorldClim tile 1...")
  if(is.null(folder)) folder <- paste0(getwd(), "/tmp")
  wclim <- geodata::worldclim_tile(var = "bio", lon = r_xy[1, 1], lat = r_xy[1, 2], path = folder)
  if(length(r_nums[[1]]) > 1){
    for(i in 2:length(r_nums[[1]])){
      message(paste0("Downloading WorldClim tile ", i, "..."))
      wc <- geodata::worldclim_tile(var = "bio", lon = r_xy[i, 1], lat = r_xy[i, 2], path = folder)
      wclim <- raster::merge(wclim, wc)
    }
  }

  # Define crop area based on buffer size
  buff_ext <- raster::extent(coords)
  ext_vals <- c()
  if(buff_ext@xmin < 0){
    ext_vals[1] <- buff_ext@xmin * (1 + buff)
  } else {
    ext_vals[1] <- buff_ext@xmin * (1 - buff)
  }
  if(buff_ext@xmax < 0){
    ext_vals[2] <- buff_ext@xmax * (1 - buff)
  } else {
    ext_vals[2] <- buff_ext@xmax * (1 + buff)
  }
  if(buff_ext@ymin < 0){
    ext_vals[3] <- buff_ext@ymin * (1 + buff)
  } else {
    ext_vals[3] <- buff_ext@ymin * (1 - buff)
  }
  if(buff_ext@ymax < 0){
    ext_vals[4] <- buff_ext@ymax * (1 - buff)
  } else {
    ext_vals[4] <- buff_ext@ymax * (1 + buff)
  }
  buff_ext <- raster::extent(ext_vals)

  # Crop raster stack to buffered area
  wclim <- raster::crop(wclim, buff_ext)

  return(wclim)
}

