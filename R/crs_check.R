
#' Check CRS of coords and layer
#'
#' @param coords sf object, data frame, or matrix representing coordinates
#' @param lyr RasterLayer or SpatRaster
#'
#' @return NULL
#'
#' @noRd
crs_check <- function(coords, lyr = NULL) {
  # Convert coords to sf
  coords <- coords_to_sf(coords)

  # Get CRS
  coords_crs <- sf::st_crs(coords)

  if (is.na(coords_crs)) warning("No CRS found for the provided coordinates. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)")

  if (!is.null(lyr)){
    lyr_crs <- sf::st_crs(lyr)
    if (is.na(lyr_crs)) warning("No CRS found for the provided raster. Make sure the coordinates and the raster have the same projection (see function details or wingen vignette)")

    if (!is.na(lyr_crs) & !is.na(coords_crs)) {
      if (coords_crs != lyr_crs) stop("CRS of the provided coordinates and raster do not match")
    }
  }
}
