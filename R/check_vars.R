
#' Check environmental layers for collinearity
#'
#' @param envlayers a RasterStack of data layers
#' @param threshold the cutoff correlation coefficient for flagging variables as collinear (numeric; defaults to 0.7)
#'
#' @return a matrix of correlation coefficients
#'
#' @export
check_env <- function(envlayers, threshold = 0.7){
  cors <- raster::layerStats(envlayers, stat = "pearson", na.rm = TRUE)
  cor_df <- cor_df_helper(cors, threshold)
  return(list(cor_df = cor_df, cor_matrix = cors))
}


#' Check extracted values for collinearity
#'
#' @param envlayers a RasterStack of data layers
#' @param coords dataframe with x and y sample coordinates
#' @param threshold the cutoff correlation coefficient for flagging variables as collinear (numeric)
#'
#' @return a matrix of correlation coefficients
#'
#' @export
check_vals <- function(envlayers, coords, threshold = 0.7){
  vals <- raster::extract(envlayers, coords)
  if(length(which(is.na(vals))) > 0) warning("NA values detected in extracted variables.")
  cors <- cor(vals, use = "na.or.complete", method = "pearson")
  cor_df <- cor_df_helper(cors, threshold)
  corrplot::corrplot.mixed(cors, upper = "ellipse", tl.pos = "lt")
  return(list(cor_df = cor_df, cor_matrix = cors))
}


#' Check geographic and environmental distances for collinearity
#'
#' @param envlayers a RasterStack of data layers
#' @param coords dataframe with x and y sample coordinates
#' @param sig significance threshold for Mantel test
#' @param type the type of geographic distance to be calculated; options are "Euclidean" for direct distance, "topographic" for topographic distances, and "resistance" for resistance distances
#' @param lyr DEM raster for calculating topographic distances or resistance raster for calculating resistance distances
#'
#' @return a list with (1) a dataframe of significantly correlated variables, (2) a matrix of p-values, (3) a matrix of Mantel's r
#'
#' @export
check_dists <- function(envlayers, coords, type = "Euclidean", lyr = NULL, sig = 0.05){
  valsNA <- raster::extract(envlayers, coords)

  if(sum(!complete.cases(valsNA))){
    warning("removing ", sum(!complete.cases(valsNA)), " locations with environmental NA values for Mantel test")
    vals <- valsNA[complete.cases(valsNA),]
    coords <- coords[complete.cases(valsNA),]
  }

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

  mantel_df <-
    data.frame(var1 = rownames(p)[row(p)[upper.tri(p)]],
               var2 = rownames(p)[col(p)[upper.tri(p)]],
               p = p[upper.tri(p)],
               r = r[upper.tri(r)]) %>%
    dplyr::filter(p < sig)

  if(nrow(mantel_df) == 1) message(paste0("Warning: The distances for 1 pair of variables are significantly correlated. algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding."))
  else if(nrow(mantel_df) > 1) message(paste0("Warning: The distances for ", nrow(mantel_df), " pairs of variables are significantly correlated. algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding."))
  return(list(mantel_df = mantel_df, p = p, r = r))
}


#' Helper function to create correlation dataframe from matrix and filter based on threshold
#'
cor_df_helper <- function(cor, threshold){
  cor_df <-
    data.frame(var1 = rownames(cors)[row(cors)[upper.tri(cors)]],
               var2 = colnames(cors)[col(cors)[upper.tri(cors)]],
               r = cors[upper.tri(cors)]) %>%
    dplyr::filter(r > threshold)

  if(nrow(cor_df) == 1) message(paste0("Warning: The extracted values for 1 pair of variables had a correlation coefficient > ", threshold, ". algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding."))
  else if(nrow(cor_df) > 1) message(paste0("Warning: The extracted values for ", nrow(cor_df), " pairs of variables had correlation coefficients > ", threshold, ". algatr recommends reducing collinearity by removing correlated variables or performing a PCA before proceeeding."))

  return(cor_df)
}

