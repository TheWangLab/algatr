#' Convert from matrix, data frame, or sf to sf (sf is a pass through)
#'
#' @param coords sf object, data frame, or matrix representing coordinates
#'
#' @return converted coords in sf format
#' @export
#'
coords_to_sf <- function(coords) {
  if (inherits(coords, "sf")) {
    return(coords)
  }
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  coords <- sf::st_as_sf(coords, coords = c("x", "y"))
  return(coords)
}

#' Convert from matrix, data frame, or sf to formatted sp
#'
#' @param coords sf object, data frame, or matrix representing coordinates
#'
#' @return converted coords in sp format
#' @export
#'
coords_to_sp <- function(coords) {
  coords <- coords_to_sf(coords)
  coords <- sf::as_Spatial(coords)
  # Needs to be x + y for kriging
  colnames(coords@coords) <- c("x", "y")
  return(coords)
}

# convert from matrix/data.frame/sf to formatted df
coords_to_df <- function(coords) {
  if (inherits(coords, "sf")) coords <- data.frame(coords_to_sp(coords)) %>% dplyr::select(-optional)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  return(coords)
}

# convert from matrix/data.frame/sf to formatted matrix
coords_to_matrix <- function(coords) {
  coords <- coords_to_df(coords)
  coords <- as.matrix(coords)
  return(coords)
}
