
#  convert from matrix/data.frame/sf to sf (sf is a pass through)
coords_to_sf <- function(coords){
  if (inherits(coords, "sf")) return(coords)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  coords <- sf::st_as_sf(coords, coords = c("x", "y"))
  return(coords)
}

#  convert from matrix/data.frame/sf to sp
coords_to_sp <- function(coords){
  coords <- coords_to_sf(coords)
  coords <- sf::as_Spatial(coords)
  # needs to be x + y for kriging
  colnames(coords@coords) <- c("x", "y")
  return(coords)
}

# convert from matrix/data.frame/sf to df
coords_to_df <- function(coords){
  if (inherits(coords, "sf")) coords <- data.frame(coords_to_sp(coords)) %>% dplyr::select(-optional)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  return(coords)
}
