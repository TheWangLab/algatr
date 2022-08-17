
#' Calculate environmental distance between environmental vars
#'
#' @param env dataframe with environmental values for each coordinate
#'
#' @return
#' @export
#'
#' @examples
env_dist <- function(env){

  # Standardize environmental variables
  scalenv <- scale(env, center = TRUE, scale = TRUE)

  distmat <- as.matrix(dist(scalenv, diag = TRUE, upper = TRUE))

  return(distmat)
}


#' Calculate geographic distance between sampling coordinates
#'
#' @param coords dataframe with x and y coordinates (MUST BE CALLED X AND Y)
#'
#' @return
#' @export
#'
#' @examples
geo_dist <- function(coords){
  # Calculate geodesic distance between points
  distmat <- geodist::geodist(coords, measure = "geodesic")

  return(distmat)
}
