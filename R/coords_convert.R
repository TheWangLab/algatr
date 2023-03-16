
#' Convert from matrix, data frame, or sf to sf
#'
#' @param coords sf object, data frame, or matrix representing coordinates
#'
#' @return converted coords in sf format
#' @export
#'
coords_to_sf <- function(coords){
  if (inherits(coords, "sf")) return(coords)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  coords <- sf::st_as_sf(coords, coords = c("x", "y"))
  return(coords)
}

#' Convert from matrix, data frame, or sf to sp
#'
#' @param coords sf object, data frame, or matrix representing coordinates
#'
#' @return converted coords in sp format
#' @export
#'
coords_to_sp <- function(coords){
  coords <- coords_to_sf(coords)
  coords <- sf::as_Spatial(coords)
  # Needs to be x + y for kriging
  colnames(coords@coords) <- c("x", "y")
  return(coords)
}
