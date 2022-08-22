

#' Calculate distance between environmental vars
#'
#' @param env dataframe with environmental values for each coordinate
#' @param stdz if TRUE then environmental values will be standardized (defaults to TRUE)
#'
#' @return
#' @export
#'
#' @examples
env_dist <- function(env, stdz = TRUE){
  if(!is.null(dim)) distmat <- dplyr::as_tibble(env) %>%  purrr::map(env_dist_helper, stdz = stdz)
  if(is.null(dim)) distmat <- env_dist_helper(env, stdz)
  return(distmat)
}

#' Helper function to convert an environemntal vector to a distance matrix
#'
#' @inheritParams env_dist
#'
#' @export
#' @noRd
env_dist_helper <- function(env, stdz = TRUE){

  # Standardize environmental variables
  if(stdz) env <- scale(env, center = TRUE, scale = TRUE)

  distmat <- as.matrix(dist(env, diag = TRUE, upper = TRUE))

  return(distmat)
}

#' Calculate geographic distance between sampling coordinates
#'
#' @param coords dataframe with x and y coordinates
#' @param type The type of geographic distance to be calculated; options are "Euclidean" for direct distance, "topographic" for topographic distances, and "resistance" for resistance distances.
#' @param lyr DEM raster for calculating topographic distances or resistance raster for calculating resistance distances
#' @details
#' Euclidean, or linear, distances are calculated using the geodist package: Padgham M, Sumner M (2021). geodist: Fast, Dependency-Free Geodesic Distance Calculations. R package version 0.0.7, Available: https://CRAN.R-project.org/package=geodist.
#' Topographic distances are calculated using the topoDistance package: Wang I.J. (2020) Topographic path analysis for modeling dispersal and functional connectivity: calculating topographic distances using the TOPODISTANCE R package. Methods in Ecology and Evolution, 11: 265-272.
#' Resistance distances are calculated using the gdistance package: van Etten, J. (2017). R package gdistance: Distances and routes on geographical grids. Journal of Statistical Software, 76(1), 1–21.
#' @return A distance matrix
#' @export
#'
#' @examples
geo_dist <- function(coords, type = "Euclidean", lyr = NULL){
  if(type == "Euclidean" | type == "euclidean" | type == "linear"){
    # Calculate geodesic distance between points
    distmat <- geodist::geodist(coords, measure = "geodesic")
  }
  else if(type == "topo" | type == "topographic"){
    if(is.null(lyr)) stop("Calculating topographic distances requires a DEM layer for argument lyr.")
    message("Calculating topo distances... This can be time consuming with many points and large rasters.")
    distmat <- topoDistance::topoDist(lyr, coords, paths = FALSE)
  }
  else if(type == "resistance" | type == "cost" | type == "res"){
    if(is.null(lyr)) stop("Calculating resistance distances requires a resistance surface for argument lyr.")
    message("Calculating resistance distances... This can be time consuming with many points and large rasters.")
    # Convert resistance surface to conductance surface
    cond.r <- 1 / lyr
    trSurface <- gdistance::transition(cond.r, transitionFunction = mean, directions = 8) # Create transition surface
    trSurface <- gdistance::geoCorrection(trSurface, type = "c", scl = FALSE)
    sp <- sp::SpatialPoints(coords = coords)
    distmat <- as.matrix(gdistance::commuteDistance(trSurface, sp)) # Calculate circuit distances
  }
  return(distmat)
}
