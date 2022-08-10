
env_dist <- function(envvals){
  # Standardize environmental variables
  scalenv <- scale(envvals, center = TRUE, scale = TRUE) 
  
  distmat <- as.matrix(dist(scalenv, diag = TRUE, upper = TRUE))
  
  return(distmat)
}


geo_dist <- function(coords){
  # calculate geodesic distance between points
  distmat <- geodist::geodist(coords, measure = "geodesic")
  
  return(distmat)
}
