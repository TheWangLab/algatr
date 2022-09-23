#' Check environmental layers for collinearity
#' @param envlayers A RasterStack of data layers
#' @param threshold The cutoff correlation coefficient for flagging variables as collinear (numeric)
#'
#' @return A matrix of correlation coefficients
#'
#' @export
check_env <- function(envlayers, threshold = 0.5){
  cors <- raster::layerStats(envlayers, stat = "pearson", na.rm = TRUE)
  cors.m <- as.matrix(cors[[1]])
  counter <- 0
  for(i in 1:(nrow(cors.m) - 1)){
    for(j in (i + 1):ncol(cors.m)){
      if(abs(cors.m[i, j]) >= threshold){
        message(paste0("Layers ", rownames(cors.m)[i], " and ", colnames(cors.m)[i], " have a correlation coefficient of ", cors.m[i, j], "."))
        counter <- counter + 1
      }
    }
  }
  if(counter == 1) message(paste0("Warning: 1 pair of layers had a correlation coefficient > ", threshold, ". algatr recommends reducing collinearity by removing correlated variables or performing raster PCA before proceeeding."))
  else if(counter > 1) message(paste0("Warning: ", counter, " pairs of layers had correlation coefficients > ", threshold, ". algatr recommends reducing collinearity by removing correlated variables or performing raster PCA before proceeeding."))
  return(cors.m)
}

#' Check extracted values for collinearity
#' @param envlayers A RasterStack of data layers
#' @param coords Dataframe with x and y sample coordinates.
#' @param threshold The cutoff correlation coefficient for flagging variables as collinear (numeric)
#'
#' @return A matrix of correlation coefficients
#'
#' @export
check_vals <- function(envlayers, coords, threshold = 0.5){
  vals <- raster::extract(envlayers, coords)
  if(length(which(is.na(vals))) > 0) warning("NA values detected in extracted variables.")
  cors <- cor(vals, use = "na.or.complete", method = "pearson")
  counter <- 0
  for(i in 1:(nrow(cors) - 1)){
    for(j in (i + 1):ncol(cors)){
      if(abs(cors[i, j]) >= threshold){
        message(paste0("Variables ", rownames(cors)[i], " and ", colnames(cors)[i], " have a correlation coefficient of ", cors[i, j], "."))
        counter <- counter + 1
      }
    }
  }
  if(counter == 1) message(paste0("Warning: The extracted values for 1 pair of variables had a correlation coefficient > ", threshold, ". algatr recommends reducing collinearity by removing correlated variables or performing PCA before proceeeding."))
  else if(counter > 1) message(paste0("Warning: The extracted values for ", counter, " pairs of variables had correlation coefficients > ", threshold, ". algatr recommends reducing collinearity by removing correlated variables or performing PCA before proceeeding."))
  return(cors)
}

#' Check geographic and environmental distances for collinearity
#' @param envlayers A RasterStack of data layers
#' @param coords Dataframe with x and y sample coordinates.
#' @param type The type of geographic distance to be calculated; options are "Euclidean" for direct distance, "topographic" for topographic distances, and "resistance" for resistance distances.
#' @param lyr DEM raster for calculating topographic distances or resistance raster for calculating resistance distances
#'
#' @return A list with matrices of p-values and Mantel's r
#'
#' @export
check_dists <- function(envlayers, coords, type = "Euclidean", lyr = NULL){
  vals <- raster::extract(envlayers, coords)
  edists <- env_dist(vals)
  gdist <- geo_dist(coords, type = type, lyr = lyr)
  gdist <- list(geo = gdist)
  all_dists <- c(gdist, edists)
  p <- r <- matrix(nrow = length(all_dists), ncol = length(all_dists))
  rownames(p) <- rownames(r) <- names(all_dists)

  for(i in 1:(length(all_dists) - 1)){
    for(j in (i+1):length(all_dists)){
      mantel_result <- vegan::mantel(all_dists[[i]], all_dists[[j]], permutations = 99)
      p[i, j] <- mantel_result$signif
      r[i, j] <- mantel_result$statistic
    }
  }

  counter <- 0
  for(i in 1:(nrow(p) - 1)){
    for(j in (i + 1):ncol(p)){
      if(p[i, j] <= 0.05){
        message(paste0(rownames(p)[i], " distances and ", rownames(p)[j], " are significantly correlated (p = ", p[i,j], " Mantel's r = ", r[i, j], "."))
        counter <- counter + 1
      }
    }
  }
  if(counter == 1) message(paste0("Warning: The distances for 1 pair of variables are sifnicantly correlated. algatr recommends reducing collinearity by removing correlated variables or performing PCA before proceeeding."))
  else if(counter > 1) message(paste0("Warning: The distances for ", counter, " pairs of variables are significantly correlated. algatr recommends reducing collinearity by removing correlated variables or performing PCA before proceeeding."))
  return(list(p, r))
}

