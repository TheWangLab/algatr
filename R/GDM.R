
#' GDM function to do everything (fit model, get coefficients, make and save raster)
#'
#' @param gendist matrix of genetic distances (must range between 0 and 1 or set scale_gendist = TRUE)
#' @param coords dataframe with x (i.e., longitude) and y (i.e., latitude) coordinates; must be in this order
#' @param envlayers envlayers for mapping (if env is provided, the dataframe column names and envlayers layer names should be the same)
#' @param env dataframe or raster object with environmental values for each coordinate; if not provided, it will be calculated based on coords/envlayers
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); (defaults to "best")
#' @param sig alpha value for significance threshold (defaults to 0.05); only used if model = "best"
#' @param nperm number of permutations to use to calculate variable importance; only used if model = "best" (defaults to 50)
#' @param geodist_type the type of geographic distance to be calculated; options are "Euclidean" (default) for direct distance, "topographic" for topographic distances, and "resistance" for resistance distances. Note: creation and plotting of the GDM raster is only possible for "Euclidean" distances
#' @param dist_lyr DEM raster for calculating topographic distances or resistance raster for calculating resistance distances
#' @param scale_gendist whether to scale genetic distance data from 0 to 1 (defaults to FALSE)
#' @param plot_vars whether to create variable vector loading plot (defaults to TRUE)
#' @param quiet whether to print output tables and figures (defaults to FALSE)
#'
#' @details
#' GDM is run using the gdm package: Fitzpatrick, M., Mokany, K., Manion, G., Nieto-Lugilde, D., & Ferrier, S. (2022). gdm: Generalized dissimilarity modeling. R package version 1.5.0-3.
#'
#' @return list with final model, predictor coefficients, and PCA RGB map
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_do_everything <- function(gendist, coords, envlayers = NULL, env = NULL, model = "best", sig = 0.05, nperm = 50,
                              geodist_type = "Euclidean", dist_lyr = NULL, scale_gendist = FALSE, plot_vars = TRUE,
                              quiet = FALSE){

  # Check CRS of envlayers and coords
  crs_check(coords, envlayers)

  # If coords not provided, make env dataframe from layers and coords
  if(is.null(env)) env <- terra::extract(envlayers, coords, ID = FALSE)

  # Run model with all defaults
  gdm_result <- gdm_run(gendist, coords = coords, env = env, model = model, sig = sig, nperm = nperm, scale_gendist = scale_gendist, geodist_type = geodist_type, dist_lyr = dist_lyr)

  # If mod is null, exit
  if(is.null(gdm_result$model)){warning("GDM model is NULL, returning NULL object"); return(NULL)}

  # Get coefficients from models and print table if specified
  coeff_df <- gdm_df(gdm_result)
  if(!quiet) print(gdm_table(gdm_result))

  # Plot I-splines if output printed
  if(!quiet) gdm_plot_isplines(gdm_result$model)

  # Create and plot map
  if(geodist_type == "Euclidean" | is.null(envlayers)) map <- gdm_map(gdm_result$model, envlayers, coords, plot_vars = plot_vars, quiet = quiet)

  # Create list to store results
  results <- list()
  # Add model
  results[["model"]] <- gdm_result$model
  # Add predictors
  results[["coeff_df"]] <- coeff_df
  # Add varimp
  results[["varimp"]] <- gdm_result$varimp
  # Add raster(s)
  if(geodist_type == "Euclidean" | is.null(envlayers)) results[["rast"]] <- map

  return(results)
}


#' Run GDM and return model object
#'
#' @inheritParams gdm_do_everything
#'
#' @family GDM functions
#'
#' @return GDM model
#' @export
#'
#' @examples
gdm_run <- function(gendist, coords, env, model = "best", sig = 0.05, nperm = 50, scale_gendist = FALSE,
                    geodist_type = "Euclidean", distPreds = NULL, dist_lyr = NULL){

  # FORMAT DATA ---------------------------------------------------------------------------------------------------

  # Extract environmental data if env is a RasterStack
  if(inherits(env, "Raster")) env <- terra::extract(env, coords, ID = FALSE)

  # Scale genetic distance data from 0 to 1
  if(scale_gendist){gendist <- scale01(gendist)}
  if(!scale_gendist & max(gendist) > 1) stop("Maximum genetic distance is greater than 1, set scale = TRUE to rescale from 0 to 1")
  if(!scale_gendist & min(gendist) < 0) stop("Minimum genetic distance is less than 0, set scale = TRUE to rescale from 0 to 1")

  # Vector of sites (for individual-based sampling, this is just assigning 1 site to each individual)
  site <- 1:nrow(gendist)

  # Bind vector of sites with gen distances
  gdmGen <- cbind(site, gendist)

  # Convert coords to df
  coords_df <- coords_to_df(coords)

  # Create dataframe of predictor variables
  gdmPred <- data.frame(site = site,
                        x = coords_df$x,
                        y = coords_df$y,
                        env)

  # Format data for GDM
  if (geodist_type == "resistance" | geodist_type == "topographic"){
    distmat <- geo_dist(coords, type = geodist_type, lyr = dist_lyr)
    gdmDist <- cbind(site, distmat)
    gdmData <- gdm::formatsitepair(gdmData, 4, predData = gdmPred, siteColumn = "site", distPreds = list(geodist = as.matrix(gdmDist)))
  } else {
    gdmData <- gdm::formatsitepair(gdmGen, bioFormat = 3, XColumn = "x", YColumn = "y", siteColumn = "site", predData = gdmPred)
  }


  # RUN GDM -------------------------------------------------------------------------------------------------------

  # If model = "full", the final GDM model is just the full model
  if(model == "full"){
    # Remove any remaining incomplete cases
    cc <- stats::complete.cases(gdmData)
    if(!all(cc)){gdmData <- gdmData[cc, ]; warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))}

    # Run GDM with all predictors
    if(geodist_type == "resistance" | geodist_type == "topographic"){
      gdm_model_final <- gdm::gdm(gdmData, geo = FALSE)
    } else {
      gdm_model_final <- gdm::gdm(gdmData, geo = TRUE)
    }
  }

  # If model = "best", conduct variable selection procedure
  if(model == "best"){
    # Add distance matrix separately
    if(geodist_type == "Euclidean"){
      geodist <- geo_dist(coords)
      gdmDist <- cbind(site, geodist)
      gdmData <- gdm::formatsitepair(gdmData, 4, predData = gdmPred, siteColumn = "site", distPreds = list(geodist = as.matrix(gdmDist)))
    }

    # Remove any remaining incomplete cases
    cc <- stats::complete.cases(gdmData)
    if(!all(cc)){gdmData <- gdmData[cc, ]; warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))}

    # Get subset of variables for final model
    gdm_varimp <- gdm_var_select(gdmData, sig = sig, nperm = nperm)
    finalvars <- gdm_varimp$finalvars

    # Stop if there are no significant final variables
    if(is.null(finalvars) | length(finalvars) == 0){
      warning("No significant combination of variables, found returning NULL object")
      return(NULL)
    }

    # Check if x is in finalvars (i.e., if geography is significant/should be included)
    if("geo" %in% finalvars){
      geo <- TRUE
      # Remove geo from finalvars before subsetting
      finalvars <- finalvars[which(finalvars != "geo")]
    } else {
      geo <- FALSE
    }

    # Subset predictor data frame
    gdmPred_final <- gdmPred[,c("site", "x", "y", finalvars)]

    # Reformat for GDM
    gdmData_final <- gdm::formatsitepair(gdmGen,
                                    bioFormat = 3,
                                    predData = gdmPred_final,
                                    XColumn = "x",
                                    YColumn = "y",
                                    siteColumn = "site")

    # Remove any remaining incomplete cases (there shouldn't be any at this point, but added as a check)
    cc <- stats::complete.cases(gdmData_final)
    if(!all(cc)){
      gdmData_final <- gdmData_final[cc, ]
      warning(paste(sum(!cc), "NA values found in final gdmData, removing;", sum(cc), "values remain"))
    }

    # Run final model
    gdm_model_final <- gdm::gdm(gdmData_final, geo = geo)

    return(list(model = gdm_model_final, pvalues = gdm_varimp$pvalues, varimp = gdm_varimp$varimp))

  }
  return(list(model = gdm_model_final, pvalues = NULL, varimp = NULL))
}



#' Get best set of variables from a GDM model
#'
#' @param gdmData data formatted using GDM package
#' @param sig sig level for determining variable significance
#' @param nperm number of permutations to run for variable testing
#'
#' @return
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_var_select <- function(gdmData, sig = 0.05, nperm = 10){
  # Check var importance/significance (THIS STEP CAN TAKE A WHILE)
  vars <- gdm::gdm.varImp(gdmData,
                     geo = FALSE,
                     splines = NULL,
                     nPerm = nperm,
                     predSelect = TRUE)

  # Get p-values from variable selection model
  pvalues <- vars[[3]]

  # Identify which cells have p-values lower than sig and not NA
  # NB: NA occurs when variables are removed during model testing
  cond <- (pvalues < sig) | is.na(pvalues)

  # Identify which columns (i.e., models) have all significant p-values (or NA)
  mods <- apply(cond, 2, all)

  # Stop if there are no models with all sig p-values
  if(all(!mods)) {warning("No significant model variable set found, returning NULL"); return(NULL)}

  # Identify the first model (i.e., minimum) that has all significant p-values
  finalmod <- min(which(mods))

  # Subset out final mod variables
  finalmod <- pvalues[,finalmod]

  # Get final variable names (i.e., names that are not NA)
  finalvars <- rownames(pvalues)[which(!is.na(finalmod))]
  finalpval <- finalmod[which(!is.na(finalmod))]

  # If the geodist matrix (matrix_1) is a significant variable, add geo to the list of vars
  if("matrix_1" %in% finalvars) {
    # Remove matrix
    finalvars <- finalvars[which(finalvars != "matrix_1")]
    # Add geo
    finalvars <- c("geo", finalvars)
  }

  return(list(finalvars = finalvars, pvalues = finalpval, varimp = vars))
}



#' Make map from model
#'
#' @param gdm_model GDM model
#' @param envlayers stack of raster layers (NAMES MUST CORRESPOND WITH GDM MODEL)
#' @param plot_vars whether to create PCA plot to help in variable and map interpretation
#' @param coords data frame with x and y coordinates
#' @param scl constant for rescaling variable vectors for plotting
#' @param display_axes display PC axes text, labels, and ticks (defaults to FALSE)
#' @inheritParams gdm_do_everything
#'
#' @return GDM RGB map
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_map <- function(gdm_model, envlayers, coords, plot_vars = TRUE, scl = 1, display_axes = FALSE, quiet = FALSE){

  # convert envlayers to SpatRaster
  if (!inherits(envlayers, "SpatRaster")) envlayers <- terra::rast(envlayers)

  # convert coords to df
  coords <- coords_to_df(coords)

  # CHECK that all of the model variables are included in the stack of environmental layers
  # Create list of environmental predictors (everything but Geographic)
  check_geo <- gdm_model$predictors == "Geographic"
  if(any(check_geo)){model_vars <- gdm_model$predictors[-which(check_geo)]} else {model_vars <- gdm_model$predictors}

  # Check that model variables are included in names of envlayers
  var_check <- model_vars %in% names(envlayers)

  # Print error with missing layers
  if(!all(var_check)){stop(paste("missing model variable(s) from raster stack:",  model_vars[!var_check]))}

  # Subset envlayers to only include variables in final model
  envlayers_sub <- terra::subset(envlayers, model_vars)


  # CREATE MAP ----------------------------------------------------------------------------------------------------

  # Transform GIS layers
  # Convert envlayers to raster
  envlayers_sub <- raster::stack(envlayers_sub)
  rastTrans <- gdm::gdm.transform(gdm_model, envlayers_sub)
  rastTrans <- terra::rast(rastTrans)

  # Remove NA values
  rastDat <- na.omit(terra::values(rastTrans))

  # Run PCA
  pcaSamp <- stats::prcomp(rastDat)

  # Count number of layers
  n_layers <- terra::nlyr(rastTrans)
  # Max number of layers to plot is 3, so adjust n_layers accordingly
  if (n_layers > 3){n_layers <- 3}

  # Make PCA raster
  pcaRast <- terra::predict(rastTrans, pcaSamp, index=1:n_layers)

  # Scale rasters to get colors (each layer will correspond with R, G, or B in the final plot)
  pcaRastRGB <- stack_to_rgb(pcaRast)

  # If there are fewer than 3 n_layers (e.g., <3 variables), the RGB plot won't work (because there isn't an R, G, and B)
  # To get around this, create a blank raster (i.e., a white raster), and add it to the stack
  if(n_layers < 3){
    warning("Fewer than three non-zero coefficients provided, adding white substitute layers to RGB plot")
    # Create white raster by multiplying a layer of pcaRast by 0 and adding 255
    white_raster <- pcaRastRGB[[1]]*0 + 255
  }

  # If n_layers = 2, you end up making a bivariate map
  if(n_layers == 2){pcaRastRGB <- c(pcaRastRGB, white_raster)}

  # If n_layers = 1, you end up making a univariate map
  if(n_layers == 1){pcaRastRGB <- c(pcaRastRGB, white_raster, white_raster)}

  # Plot raster if quiet = FALSE
  if(!quiet) terra::plotRGB(pcaRastRGB, r = 1, g = 2, b = 3)
  if(!is.null(coords)) points(coords, cex = 1.5)

  # Plot variable vectors
  if(plot_vars & (n_layers == 3)){
    # TODO [EAC] need to add quiet here?
    gdm_plot_vars(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = scl, display_axes = display_axes, quiet = quiet)
  }

  if(plot_vars & (n_layers != 3)){
    warning("variable vector plot is not available for model with fewer than 3 final variables, skipping...")
  }

  s <- list(rastTrans, pcaRastRGB)
  names(s) <- c("rastTrans", "pcaRastRGB")
  return(s)
}


#' Plot I-splines for each variable
#'
#' @param gdm_model GDM model
#'
#' @return plot for each I-spline
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_plot_isplines <- function(gdm_model){
  gdm_model_splineDat <- gdm::isplineExtract(gdm_model)

  purrr::walk(1:ncol(gdm_model_splineDat$x), function(i){
    dat <- cbind(as.data.frame(gdm_model_splineDat$x[,i]), as.data.frame(gdm_model_splineDat$y[,i]))
    plot <- ggplot2::ggplot(dat) +
      ggplot2::geom_line(ggplot2::aes(x = gdm_model_splineDat$x[,i], y = gdm_model_splineDat$y[,i])) +
      ggplot2::theme_bw() +
      ggplot2::xlab(colnames(gdm_model_splineDat$x)[i]) +
      ggplot2::ylab("Partial Regression Distance")

    print(plot)})
  }



#' Plot compositional dissimilarity spline plots
#'
#' @description generates two plots: a plot of the observed response data against raw ecological distance from the model, and a plot of the observed response against the predicted response from the model (after link function is applied)
#' @param gdm_model GDM model
#'
#' @return two spline plots of compositional dissimilarity
#'
#' @family GDM functions
#' @details code is modified from the `plot.gdm()` function in the gdm package (Fitzpatrick et al. 2022)
#'
#' @export
#'
#' @examples
gdm_plot_diss <- function(gdm_model){
  obs <- tidyr::as_tibble(gdm_model$observed) %>% dplyr::rename(observed = value)
  pred <- tidyr::as_tibble(gdm_model$predicted) %>% dplyr::rename(predicted = value)
  ecol <- tidyr::as_tibble(gdm_model$ecological) %>% dplyr::rename(ecological = value)

  dat <- cbind(obs, pred, ecol)
  datL <- nrow(dat)

  # Get data for overlaid lines
  overlayX_ecol <- seq(from = min(dat$ecological), to = max(dat$ecological), length = datL)
  overlayY_ecol <- 1-exp(-overlayX_ecol)
  overlayX_pred <- seq(from = min(dat$predicted), to = max(dat$predicted), length = datL)
  overlayY_pred <- 1-exp(-overlayX_pred)

  plot_ecol <-
    ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(x = ecological, y = observed), color = "darkgrey", alpha = 0.6, size = 2) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    ggplot2::xlab("Predicted ecological distance") +
    ggplot2::ylab("Observed compositional dissimilarity") +
    ggplot2::geom_line(ggplot2::aes(x = overlayX_ecol, y = overlayY_ecol), size = 1)

  plot_pred <-
    ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(x = predicted, y = observed), color = "darkgrey", alpha = 0.6, size = 2) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    ggplot2::xlab("Predicted compositional dissimilarity") +
    ggplot2::ylab("Observed compositional dissimilarity") +
    ggplot2::geom_line(ggplot2::aes(x = overlayX_pred, y = overlayY_pred), size = 1)

  plot <- cowplot::plot_grid(plot_ecol, plot_pred, nrow = 1)

  print(plot)
}



#' Create a PCA plot for GDM
#'
#' @param pcaSamp PCA results from running prcomp()
#' @param pcaRast raster PCA
#' @param pcaRastRGB raster PCA rescaled to RGB
#' @param coords dataframe with x and y coordinates
#' @param x x-axis PC
#' @param y y-axis PC
#' @param scl constant for rescaling variable vectors for plotting
#' @param display_axes whether to display axes
#' @inheritParams gdm_do_everything
#'
#' @return GDM PCA plot
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_plot_vars <- function(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = 1, display_axes = FALSE, quiet = FALSE){

  # Confirm there are exactly 3 axes
  if(terra::nlyr(pcaRastRGB) > 3){stop("Only three PC layers (RGB) can be used for creating the variable plot (too many provided)")}
  if(terra::nlyr(pcaRastRGB) < 3){stop("Need exactly three PC layers (RGB) for creating the variable plot (too few provided)")}

  # GET PCA DATA ----------------------------------------------------------------------------------------------------

  # Make data frame from PC results
  xpc <- data.frame(pcaSamp$x[,1:3])

  # Get variable rotations
  varpc <- data.frame(varnames = rownames(pcaSamp$rotation), pcaSamp$rotation)

  # Get PC values for each coord
  pcavals <- data.frame(terra::extract(pcaRast, coords, ID = FALSE))
  colnames(pcavals) <- colnames(xpc)

  # Rescale var loadings with individual loadings so they fit in the plot nicely
  scldat <- min(
    (max(pcavals[,y], na.rm = TRUE) - min(pcavals[,y], na.rm = TRUE)/(max(varpc[,y], na.rm = TRUE)-min(varpc[,y], na.rm = TRUE))),
    (max(pcavals[,x], na.rm = TRUE) - min(pcavals[,x], na.rm = TRUE)/(max(varpc[,x], na.rm = TRUE)-min(varpc[,x], na.rm = TRUE)))
  )

  # Additionally use a constant scale val (scl) to shrink the final vectors (again for plotting nicely)
  varpc <- data.frame(varpc, v1 = scl * scldat * varpc[,x],
                     v2 = scl * scldat *  varpc[,y])


  # GET RGB VALS FOR EACH COORD----------------------------------------------------------------------------------------

  pcavalsRGB <- data.frame(terra::extract(pcaRastRGB, coords, ID = FALSE))
  colnames(pcavalsRGB) <- colnames(xpc)

  # Create vector of RGB colors for plotting
  pcacols <- apply(pcavalsRGB, 1, create_rgb_vec)

  # GET RGB VALS FOR ENTIRE RASTER-------------------------------------------------------------------------------------

  # Get sample
  s <- sample(1:terra::ncell(pcaRast), 10000)

  # Get all PC values from raster and remove NAs
  rastvals <- data.frame(values(pcaRast))[s,]
  colnames(rastvals) <- colnames(xpc)
  rastvals <- rastvals[stats::complete.cases(rastvals),]

  # Get all RGB values from raster and remove NAs
  rastvalsRGB <- data.frame(values(pcaRastRGB))[s,]
  colnames(rastvalsRGB) <- colnames(rastvals)
  rastvalsRGB <- rastvalsRGB[stats::complete.cases(rastvalsRGB),]

  # Create vector of RGB colors for plotting
  rastpcacols <- apply(rastvalsRGB, 1, create_rgb_vec)

  # FINAL PLOT----------------------------------------------------------------------------------------------------------

  # Build base plot
  # Plot points colored by RGB with variable vectors
  plot <- ggplot2::ggplot() +

    # Create axes that cross through origin
    {if(display_axes)ggplot2::geom_hline(yintercept = 0, size=0.2, col = "gray")} +
    {if(display_axes)ggplot2::geom_vline(xintercept = 0, size=0.2, col = "gray")} +

    # Plot points from entire raster
    ggplot2::geom_point(data = rastvals, ggplot2::aes_string(x = x, y = y), col = rastpcacols, size = 4, alpha = 0.02) +

    # Plot coord values
    ggplot2::geom_point(data = pcavals, ggplot2::aes_string(x = x, y = y), fill = pcacols, col = "black", pch = 21, size = 3) +

    # Plot variable vectors
    ggplot2::geom_text(data = varpc, ggplot2::aes(x = v1, y = v2, label = varnames), size = 4, vjust = 1) +
    ggplot2::geom_segment(data = varpc, ggplot2::aes(x = 0, y = 0, xend = v1, yend = v2), arrow = ggplot2::arrow(length = ggplot2::unit(0.2,"cm"))) +

    # Plot formatting
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   aspect.ratio = 1)

  # Build plot without PC axes displayed
  if(display_axes == FALSE){
    plot <- plot +
      # Remove axes
      ggplot2::theme(axis.title = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank())
  }

  # Plot
  if(!quiet) print(plot)
}

#' Helper function to create rgb vector
#'
#' @export
#' @noRd
create_rgb_vec <- function(vec){
  if(any(is.na(vec))) x <- NA else x <- rgb(vec[1], vec[2], vec[3], maxColorValue = 255)
  return(x)
}

#' Scale a raster stack from 0 to 255
#'
#' @param s RasterStack
#'
#' @noRd
#' @export
stack_to_rgb <- function(s){
  stack_list <- as.list(s)
  new_stack <- terra::rast(purrr::map(stack_list, raster_to_rgb))
  return(new_stack)
}

#' Scale raster from 0 to 255
#'
#' @param r RasterLayer
#'
#' @noRd
#' @export
raster_to_rgb <- function(r){
  rmax <- terra::minmax(r)["max",]
  rmin <- terra::minmax(r)["min",]
  if(rmin == 0){r <- 255} else {r <- (r - rmin) / (rmax - rmin) * 255}
}


#' Get coefficients for each predictor
#'
#' @param gdm_model model of type `gdm`
#'
#' @return dataframe of coefficients for GDM
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_coeffs <- function(gdm_model){
  # Vector to store coefficient sums
  coefSums <- c()

  # Sum coefficients for each predictor (each has 3 splines which you sum)
  for (i in 1:length(gdm_model$predictors)){

    # Create starting index (* 3 b/c there are 3 splines and - 2 to get the first index)
    j <- (i * 3) - 2

    # Subset out three splines and sum them (j = first index, j + 2 = last index)
    coefSums[i] <- sum(gdm_model$coefficients[j:(j+2)])

  }

  # Add those values to a simple dataframe
  coeffs <- data.frame(predictor = gdm_model$predictors, coefficient = coefSums)

  return(coeffs)
}

#' Create dataframe of GDM results
#'
#' @param gdm_result output of \link[algatr]{gdm_run}
#'
#' @return dataframe of gdm model coefficients
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_df <- function(gdm_result){
  coeff_df <- gdm_coeffs(gdm_result$model)
  if(!is.null(gdm_result$pvalues)) coeff_df$p <- gdm_result$pvalues
  return(coeff_df)
}

#' Create `gt` table of GDM results
#'
#' @param gdm_result output of \link[algatr]{gdm_run} or \link[algatr]{gdm_do_everything}
#' @param digits number of digits to include (defaults to 2)
#'
#' @family GDM functions
#'
#' @return An object of class `gt_tbl`
#' @export
gdm_table <- function(gdm_result, digits = 2, summary_stats = TRUE, footnote = TRUE){

  gdm_df <- gdm_df(gdm_result)

  d <- max(abs(min(gdm_df$coefficient, na.rm = TRUE)), abs(max(gdm_df$coefficient)))

  suppressWarnings({
    tbl <- gdm_df  %>%
      gt::gt() %>%
      gtExtras::gt_hulk_col_numeric(coefficient, trim = TRUE, domain = c(-d,d)) %>%
      gt::sub_missing(missing_text = "")

    if (summary_stats) {

      stat_names <- c("% Explained:")
      stats <- as.numeric(gdm_result$model$explained)
      gdm_df <- gdm_df %>%
        rbind(purrr::map2_dfr(.x = stat_names, .y = stats, .f = make_stat_vec, gdm_df)) %>%
        dplyr::mutate(dplyr::across(-c(predictor), as.numeric))


      tbl <- gdm_df %>%
        gt::gt() %>%
        gtExtras::gt_hulk_col_numeric(coefficient, trim = TRUE, domain = c(-d,d)) %>%
        gt::sub_missing(missing_text = "") %>%
        gt::tab_row_group(label = NA, id = "model", rows = which(!(gdm_df$predictor %in% stat_names))) %>%
        gtExtras::gt_highlight_rows(rows = which(gdm_df$predictor %in% stat_names), fill = "white") %>%
        gt::tab_style(
          style = list(gt::cell_borders(sides = "top", color = "white"),
                       gt::cell_text(align = "left")),
          locations = gt::cells_body(rows = which(gdm_df$var %in% stat_names))
        )
    }

    if (footnote & summary_stats) tbl <- tbl %>% gt::tab_footnote(footnote = "The percentage of null deviance explained by the fitted GDM model.",
                                                                  locations = gt::cells_body(rows = which(gdm_df$predictor == "% Explained:"), columns = coefficient),
                                                                  placement = "right")


    if(!is.null(digits)) tbl <- tbl %>% gt::fmt_number(columns = -predictor, decimals = 2)

  })

  tbl
}

#' Scale genetic distances from 0 to 1
#'
#' @param x genetic distance matrix
#'
#' @return genetic distance matrix scaled from 0 to 1
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
scale01 <- function(x){(x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}

