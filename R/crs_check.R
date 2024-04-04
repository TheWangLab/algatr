#' Check CRS of coords and layer
#'
#' @param coords sf object, data frame, or matrix representing coordinates
#' @param lyr SpatRaster or Raster* object
#'
#' @return NULL
#' @export
#'
#' @noRd
crs_check <- function(coords = NULL, lyr = NULL) {
  if (!is.null(coords)) {
    # convert coords to sf
    coords <- coords_to_sf(coords)
    # get CRS
    coords_crs <- sf::st_crs(coords)
    if (is.na(coords_crs)) warning("No CRS found for the provided coordinates. Make sure the coordinates and the raster have the same projection (see function details or vignette)")
  }

  if (!is.null(lyr)) {
    lyr_crs <- sf::st_crs(lyr)
    if (is.na(lyr_crs)) warning("No CRS found for the provided raster. Make sure the coordinates and the raster have the same projection (see function details or vignette)")
  }

  if (!is.null(lyr) & !is.null(coords)) {
    if (!is.na(lyr_crs) & !is.na(coords_crs)) {
      if (coords_crs != lyr_crs) stop("CRS of the provided coordinates and raster do not match")
    }
  }
}
