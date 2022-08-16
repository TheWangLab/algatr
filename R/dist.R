
env_dist <- function(env){
  # Standardize environmental variables
  scalenv <- scale(env, center = TRUE, scale = TRUE)

  distmat <- as.matrix(dist(scalenv, diag = TRUE, upper = TRUE))

  return(distmat)
}


geo_dist <- function(coords){
  # calculate geodesic distance between points
  distmat <- geodist::geodist(coords, measure = "geodesic")

  return(distmat)
}
