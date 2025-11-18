#' GDM function to do everything (fit model, get coefficients, make and save raster)
#'
#' @param gendist matrix of genetic distances (must range between 0 and 1 or set scale_gendist = TRUE)
#' @param coords dataframe with x (i.e., longitude) and y (i.e., latitude) coordinates; must be in this order
#' @param envlayers SpatRaster or Raster* object for mapping (if `env`` is provided, the dataframe column names and `envlayers`` layer names should be the same)
#' @param env dataframe or raster object with environmental values for each coordinate; if not provided, it will be calculated based on coords/envlayers
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "full"
#' @param sig alpha value for significance threshold (defaults to 0.05); only used if model = "best"
#' @param nperm number of permutations to use to calculate variable importance; only used if model = "best" (defaults to 50)
#' @param geodist_type the type of geographic distance to be calculated; options are "Euclidean" (default) for direct distance, "topographic" for topographic distances, and "resistance" for resistance distances. Note: creation and plotting of the GDM raster is only possible for "Euclidean" distances
#' @param dist_lyr DEM raster for calculating topographic distances or resistance raster for calculating resistance distances
#' @param scale_gendist whether to scale genetic distance data from 0 to 1 (defaults to FALSE)
#' @param plot_vars whether to create PCA plot to help in variable and map interpretation (defaults to TRUE)
#' @param quiet whether to operate quietly and suppress the output of tables and figures (defaults to FALSE)
#'
#' @details
#' GDM is run using the gdm package: Fitzpatrick, M., Mokany, K., Manion, G., Nieto-Lugilde, D., & Ferrier, S. (2022). gdm: Generalized dissimilarity modeling. R package version 1.5.0-3.
#'
#' @return list with final model, predictor coefficients, and PCA RGB map
#'
#' @family GDM functions
#'
#' @export
gdm_do_everything <- function(gendist, coords, envlayers = NULL, env = NULL, model = "full", sig = 0.05, nperm = 50,
                              geodist_type = "Euclidean", dist_lyr = NULL, scale_gendist = FALSE, plot_vars = TRUE,
                              quiet = FALSE) {
  message("Please be aware: the do_everything functions are meant to be exploratory. We do not recommend their use for final analyses unless certain they are properly parameterized.")

  # Check CRS of envlayers and coords
  crs_check(coords, envlayers)

  # If coords not provided, make env dataframe from layers and coords
  if (is.null(env)) env <- terra::extract(envlayers, coords, ID = FALSE)

  # Run model
  gdm_run_safely <- purrr::safely(gdm_run, quiet = FALSE)
  gdm_result <- gdm_run_safely(gendist, coords = coords, env = env, model = model, sig = sig, nperm = nperm, scale_gendist = scale_gendist, geodist_type = geodist_type, dist_lyr = dist_lyr)

  if (is.null(gdm_result$result) & model == "best") {
    warning("failed to fit best model, rerunning GDM with full model")
    gdm_result <- gdm_run_safely(gendist, coords = coords, env = env, model = "full", sig = sig, nperm = nperm, scale_gendist = scale_gendist, geodist_type = geodist_type, dist_lyr = dist_lyr)
  }

  gdm_result <- gdm_result$result

  # If mod is null, exit
  if (is.null(gdm_result$model)) {
    warning("GDM model is NULL, returning NULL object")
    return(NULL)
  }

  # Get coefficients from models and print table if specified
  coeff_df <- gdm_df(gdm_result)
  if (!quiet) print(gdm_table(gdm_result))

  # Plot I-splines if output printed
  if (!quiet) print(gdm_plot_isplines(gdm_result$model))

  # check if all env splines are zero
  zero_env <-
    coeff_df %>%
    dplyr::filter(predictor != "Geographic") %>%
    # used instead of summarize for cases where Geographic is the only variable
    dplyr::reframe(sum = coefficient) %>%
    dplyr::pull() %>%
    sum()

  # if Geographic is the only variable zero_env will be an empty vector
  # replace with 0 for the logical test
  if (length(zero_env) == 0) zero_env <- 0
  if (zero_env == 0) map <- NULL

  # Create and plot map
  if (geodist_type == "Euclidean" & !is.null(envlayers)) {
    if (zero_env == 0){
      warning("All model splines for environmental variables are zero, skipping creation of GDM map")
    } else {
      map <- gdm_map(gdm_result$model, envlayers, coords, plot_vars = plot_vars, quiet = quiet)
    }
  }

  # Create list to store results
  results <- list()
  # Add model
  results[["model"]] <- gdm_result$model
  # Add predictors
  results[["coeff_df"]] <- coeff_df
  # Add varimp
  results[["varimp"]] <- gdm_result$varimp
  # Add raster(s)
  if (geodist_type == "Euclidean" & !is.null(envlayers)) results[["rast"]] <- map else results[["rast"]] <- NULL

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
gdm_run <- function(gendist, coords, env, model = "full", sig = 0.05, nperm = 50, scale_gendist = FALSE,
                    geodist_type = "Euclidean", distPreds = NULL, dist_lyr = NULL) {
  # FORMAT DATA ---------------------------------------------------------------------------------------------------
  
  # Create GDM formatted data objects
  formatted_data <- 
    gdm_format(
      gendist = gendist, 
      coords = coords, 
      env = env,
      scale_gendist = scale_gendist, 
      geodist_type = geodist_type, 
      distPreds = distPreds, 
      dist_lyr = dist_lyr,
      gdmPred = TRUE,
      gdmGen = TRUE
      )
  
  gdmData <- formatted_data$gdmData
  gdmPred <- formatted_data$gdmPred
  gdmGen <- formatted_data$gdmGen

  # Vector of sites (for individual-based sampling, this is just assigning 1 site to each individual)
  site <- 1:nrow(gendist)

  # RUN GDM -------------------------------------------------------------------------------------------------------

  # If model = "full", the final GDM model is just the full model
  if (model == "full") {
    # Remove any remaining incomplete cases
    cc <- stats::complete.cases(gdmData)
    if (!all(cc)) {
      gdmData <- gdmData[cc, ]
      warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))
    }

    # Run GDM with all predictors
    if (geodist_type == "resistance" | geodist_type == "topographic") {
      gdm_model_final <- gdm::gdm(gdmData, geo = FALSE)
    } else {
      gdm_model_final <- gdm::gdm(gdmData, geo = TRUE)
    }
  }

  # If model = "best", conduct variable selection procedure
  if (model == "best") {
    # Add Euclidean distance matrix separately (otherwise geography will not be evaluated for removal)
    if (geodist_type == "Euclidean") {
      geodist <- geo_dist(coords)
      gdmDist <- cbind(site, geodist)
      gdmData <- gdm::formatsitepair(gdmData, 4, predData = gdmPred, siteColumn = "site", distPreds = list(geodist = as.matrix(gdmDist)))
    }

    # Remove any remaining incomplete cases
    cc <- stats::complete.cases(gdmData)
    if (!all(cc)) {
      gdmData <- gdmData[cc, ]
      warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))
    }

    # Get subset of variables for final model
    gdm_varimp <- gdm_var_sel(gdmData, sig = sig, nperm = nperm)
    finalvars <- gdm_varimp$finalvars
    # For code testing you can set finalvars to a vector (i.e. c("geo", "CA_rPCA1"))

    # Stop if there are no significant final variables
    if (is.null(finalvars) | length(finalvars) == 0) {
      warning("No significant combination of variables, found returning NULL object")
      return(NULL)
    }

    # Check if geo is in finalvars (i.e., if geography is significant/should be included)
    # geo is the name given to any geographic distance matrix (matrix_1) by the gdm_var_sel function
    if ("geo" %in% finalvars) {
      geo <- TRUE
      # Remove geo from finalvars before subsetting
      finalvars <- finalvars[which(finalvars != "geo")]
    } else {
      geo <- FALSE
    }

    # Subset predictor data frame
    gdmPred_final <- gdmPred[, c("site", "x", "y", finalvars)]
    gdmData_final <- gdm::formatsitepair(gdmGen, bioFormat = 3, XColumn = "x", YColumn = "y", siteColumn = "site", predData = gdmPred_final)

    # If geo = TRUE and geodist_type is resistance or topographic, add the gdmDist matrix 
    if (geodist_type == "resistance" | geodist_type == "topographic"){
      if (geo){
        gdmData_final <- gdm::formatsitepair(gdmData_final, 4, predData = gdmPred_final, siteColumn = "site", distPreds = list(geodist = as.matrix(gdmDist)))
        # Set geo to FALSE (since it's already included, we don't want to also include Euclidean distances)
        geo <- FALSE
      } 
    } 
   
    # Remove any remaining incomplete cases (there shouldn't be any at this point, but added as a check)
    cc <- stats::complete.cases(gdmData_final)
    if (!all(cc)) {
      gdmData_final <- gdmData_final[cc, ]
      warning(paste(sum(!cc), "NA values found in final gdmData, removing;", sum(cc), "values remain"))
    }

    # Run final model
    gdm_model_final <- gdm::gdm(gdmData_final, geo = geo)

    return(list(model = gdm_model_final, pvalues = gdm_varimp$pvalues, varimp = gdm_varimp$varimp))
  }
  return(list(model = gdm_model_final, pvalues = NULL, varimp = NULL))
}

#' Format Data for Generalized Dissimilarity Modeling (GDM)
#'
#' @inheritParams gdm_do_everything
#' @param gdmGen whether to include the gdm formatted genetic data seperately (defaults to FALSE). This dataframe contains the genetic distance matrix with an additional column for site number.
#' @param gdmPred whether to include the gdm formatted predictor data seperately (defaults to FALSE). This dataframe contains the site number, coordinates, and environmental values at each site
#'
#' @family GDM functions
#'
#' @return either a gdmData object if gdmGen and gdmPred are FALSE or a list of gdm data objects
#' @export
gdm_format <- function(gendist, coords, env, scale_gendist = FALSE, geodist_type = "Euclidean", distPreds = NULL, dist_lyr = NULL, gdmPred = FALSE, gdmGen = FALSE) {

  # Rename gdmPred/gdmGen arguments as they will later be overwritten
  gdmPred_lgl <- gdmPred
  gdmGen_lgl <- gdmGen
  
  # Convert env to spat raster if it is a RasterLayer/RasterStack
  if (inherits(env, "Raster")) env <- terra::rast(env)

  # Extract environmental data if env is a raster
  if (inherits(env, "SpatRaster")) env <- terra::extract(env, coords, ID = FALSE)

  # Scale genetic distance data from 0 to 1
  if (scale_gendist) {
    gendist <- scale01(gendist)
  }
  if (!scale_gendist & max(gendist) > 1) stop("Maximum genetic distance is greater than 1, set scale_gendist = TRUE to rescale from 0 to 1")
  if (!scale_gendist & min(gendist) < 0) stop("Minimum genetic distance is less than 0, set scale_gendist = TRUE to rescale from 0 to 1")

  # Vector of sites (for individual-based sampling, this is just assigning 1 site to each individual)
  site <- 1:nrow(gendist)

  # Bind vector of sites with gen distances
  gdmGen <- cbind(site, gendist)

  # Convert coords to df
  coords_df <- coords_to_df(coords)

  # Create dataframe of predictor variables
  gdmPred <- data.frame(
    site = site,
    x = coords_df$x,
    y = coords_df$y,
    env
  )

  # Format data for GDM
  gdmData <- gdm::formatsitepair(gdmGen, bioFormat = 3, XColumn = "x", YColumn = "y", siteColumn = "site", predData = gdmPred)

  if (geodist_type == "resistance" | geodist_type == "topographic") {
    distmat <- geo_dist(coords, type = geodist_type, lyr = dist_lyr)
    gdmDist <- cbind(site, distmat)
    gdmData <- gdm::formatsitepair(gdmData, 4, predData = gdmPred, siteColumn = "site", distPreds = list(geodist = as.matrix(gdmDist)))
  } 

 
  if (!gdmPred_lgl & !gdmGen_lgl) return(gdmData)
  if (!gdmPred_lgl) gdmPred <- NULL
  if (!gdmGen_lgl) gdmGen <- NULL
  result <- purrr::compact(list(gdmData = gdmData, gdmGen = gdmGen, gdmPred = gdmPred))

  return(result)
}

#' Get best set of variables from a GDM model
#'
#' @param gdmData data formatted using GDM package
#' @param sig sig level for determining variable significance
#' @param nperm number of permutations to run for variable testing
#'
#' @family GDM functions
#'
#' @export
gdm_var_sel <- function(gdmData, sig = 0.05, nperm = 10) {
  # Check var importance/significance (THIS STEP CAN TAKE A WHILE)
  vars <- gdm::gdm.varImp(gdmData,
    geo = FALSE,
    splines = NULL,
    nPerm = nperm,
    predSelect = TRUE
  )

  # Get p-values from variable selection model
  pvalues <- vars[[3]]

  # Extract variable names from selected model
  finalvars <- names(pvalues)

  # If the geodist matrix (matrix_1) is a significant variable, add geo to the list of vars
  if ("matrix_1" %in% finalvars) {
    # Remove matrix
    finalvars <- finalvars[which(finalvars != "matrix_1")]
    # Add geo
    finalvars <- c("geo", finalvars)
  }

  return(list(finalvars = finalvars, pvalues = pvalues, varimp = vars))
}

#' Make map from model
#'
#' @param gdm_model GDM model
#' @param envlayers SpatRaster or Raster* object (LAYER NAMES MUST CORRESPOND WITH GDM MODEL)
#' @param coords data frame with x and y coordinates
#' @param scl constant for rescaling variable vectors for plotting (defaults to 1)
#' @param display_axes display PC axes text, labels, and ticks (defaults to FALSE)
#' @inheritParams gdm_do_everything
#'
#' @return GDM RGB map
#'
#' @family GDM functions
#'
#' @export
gdm_map <- function(gdm_model, envlayers, coords, plot_vars = TRUE, scl = 1, display_axes = FALSE, quiet = FALSE) {
  # convert envlayers to SpatRaster
  if (!inherits(envlayers, "SpatRaster")) envlayers <- terra::rast(envlayers)

  # convert coords to df
  coords <- coords_to_df(coords)

  # CHECK that all of the model variables are included in the stack of environmental layers
  # Create list of environmental predictors (everything but Geographic)
  check_geo <- gdm_model$predictors == "Geographic"
  if (any(check_geo)) {
    model_vars <- gdm_model$predictors[-which(check_geo)]
  } else {
    model_vars <- gdm_model$predictors
  }

  # Check that model variables are included in names of envlayers
  var_check <- model_vars %in% names(envlayers)

  # Print error with missing layers
  if (!all(var_check)) {
    stop(paste("missing model variable(s) from raster stack:", model_vars[!var_check]))
  }

  # Subset envlayers to only include variables in final model
  envlayers_sub <- terra::subset(envlayers, model_vars)

  # CREATE MAP ----------------------------------------------------------------------------------------------------

  # Transform environmental layers
  # TEMPORARY: In new versions of GDM, the input/output rasters are SpatRasters, but for the old version they are rasters
  if (packageVersion("gdm") >= "1.6.0-4") {
    rastTrans <- gdm::gdm.transform(gdm_model, envlayers_sub)
  } else {
    envlayers_sub_raster <- raster::stack(envlayers_sub)
    rastTrans <- gdm::gdm.transform(gdm_model, envlayers_sub_raster)
    rastTrans <- terra::rast(rastTrans)
  }

  # Remove NA values
  rastDat <- na.omit(terra::values(rastTrans))

  # Run PCA
  pcaSamp <- stats::prcomp(rastDat)

  # Count number of layers
  n_layers <- terra::nlyr(rastTrans)
  # Max number of layers to plot is 3, so adjust n_layers accordingly
  if (n_layers > 3) {
    n_layers <- 3
  }

  # Check if there are only coordinate layers
  # If there are only coordinate layers (i.e., no env layers) than you cannot create the map
  if (all(names(rastTrans) %in% c("xCoord", "yCoord"))) stop("All model splines for environmental variables are zero")

  # Make PCA raster
  pcaRast <- terra::predict(rastTrans, pcaSamp, index = 1:n_layers)

  # Scale rasters to get colors (each layer will correspond with R, G, or B in the final plot)
  pcaRastRGB <- stack_to_rgb(pcaRast)

  # If there are fewer than 3 n_layers (e.g., <3 variables), the RGB plot won't work (because there isn't an R, G, and B)
  # To get around this, create a blank raster (i.e., a white raster), and add it to the stack
  if (n_layers < 3) {
    warning("Fewer than three non-zero coefficients provided, adding white substitute layers to RGB plot")
    # Create white raster by multiplying a layer of pcaRast by 0 and adding 255
    white_raster <- pcaRastRGB[[1]] * 0 + 255
  }

  # If n_layers = 2, you end up making a bivariate map
  if (n_layers == 2) {
    pcaRastRGB <- c(pcaRastRGB, white_raster)
  }

  # If n_layers = 1, you end up making a univariate map
  if (n_layers == 1) {
    pcaRastRGB <- c(pcaRastRGB, white_raster, white_raster)
  }

  # Plot raster if quiet = FALSE
  if (!quiet) terra::plotRGB(pcaRastRGB, r = 1, g = 2, b = 3)
  if (!is.null(coords) & !quiet) points(coords, cex = 1.5)

  # Plot variable vectors
  if (!quiet & plot_vars & (n_layers == 3)) {
    gdm_plot_vars(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = scl, display_axes = display_axes)
  }

  if (!quiet & plot_vars & (n_layers != 3)) {
    warning("variable vector plot is not available for model with fewer than 3 final variables, skipping...")
  }

  s <- list(rastTrans, pcaRastRGB)
  names(s) <- c("rastTrans", "pcaRastRGB")
  return(s)
}


#' Plot I-splines for each variable
#'
#' @param gdm_model GDM model
#' @param scales Whether scales should be free ("free"; default), free in one dimension ("free_x", "free_y") or fixed ("fixed"). We recommend setting this to "free_x" to allow the x-axis to vary while keeping the y-axis fixed across all plots such that relative importance can be visualized.
#' @param nrow Number of rows 
#' @param ncol Number of cols
#' 
#' @return plot for each I-spline
#'
#' @family GDM functions
#'
#' @export
gdm_plot_isplines <- function(gdm_model, scales = "free", nrow = NULL, ncol = NULL) {
  gdm_model_splineDat <- gdm::isplineExtract(gdm_model)
  
  gdm_spline_df <- 
    dplyr::bind_rows(
      data.frame(gdm_model_splineDat$x, var = "x", ID = 1:nrow(gdm_model_splineDat$x)), 
      data.frame(gdm_model_splineDat$y, var = "y", ID = 1:nrow(gdm_model_splineDat$y))
      ) %>%
    tidyr::pivot_longer(-c("var", "ID"), names_to = "name", values_to = "value") %>%
    tidyr::pivot_wider(names_from = "var", values_from = "value") %>%
    dplyr::mutate(name = factor(name, levels = unique(name)))
  
  plt <-
    ggplot2::ggplot(gdm_spline_df) +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y)) +
      ggplot2::facet_wrap(~name, scales = scales, nrow = nrow, ncol = ncol, strip.position = "bottom") +
      ggplot2::theme_bw() +
      ggplot2::xlab("") +
      ggplot2::ylab("Partial Regression Distance") +
      ggplot2::theme(strip.placement = "outside", strip.background = ggplot2::element_blank())

  return(plt)
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
gdm_plot_diss <- function(gdm_model) {
  obs <- tidyr::as_tibble(gdm_model$observed) %>% dplyr::rename(observed = value)
  pred <- tidyr::as_tibble(gdm_model$predicted) %>% dplyr::rename(predicted = value)
  ecol <- tidyr::as_tibble(gdm_model$ecological) %>% dplyr::rename(ecological = value)

  dat <- cbind(obs, pred, ecol)
  datL <- nrow(dat)

  # Get data for overlaid lines
  overlayX_ecol <- seq(from = min(dat$ecological), to = max(dat$ecological), length = datL)
  overlayY_ecol <- 1 - exp(-overlayX_ecol)
  overlayY_pred <- overlayX_pred <- seq(from = min(dat$predicted), to = max(dat$predicted), length = datL)

  plot_ecol <-
    ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(x = ecological, y = observed), color = "darkgrey", alpha = 0.6, size = 2) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::xlab("Predicted ecological distance") +
    ggplot2::ylab("Observed compositional dissimilarity") +
    ggplot2::geom_line(ggplot2::aes(x = overlayX_ecol, y = overlayY_ecol), size = 1)

  plot_pred <-
    ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(x = predicted, y = observed), color = "darkgrey", alpha = 0.6, size = 2) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
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
#'
#' @return GDM PCA plot
#'
#' @family GDM functions
#'
#' @export
gdm_plot_vars <- function(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = 1, display_axes = FALSE) {
  # Confirm there are exactly 3 axes
  if (terra::nlyr(pcaRastRGB) > 3) {
    stop("Only three PC layers (RGB) can be used for creating the variable plot (too many provided)")
  }
  if (terra::nlyr(pcaRastRGB) < 3) {
    stop("Need exactly three PC layers (RGB) for creating the variable plot (too few provided)")
  }

  # GET PCA DATA ----------------------------------------------------------------------------------------------------

  # Make data frame from PC results
  xpc <- data.frame(pcaSamp$x[, 1:3])

  # Get variable rotations
  varpc <- data.frame(varnames = rownames(pcaSamp$rotation), pcaSamp$rotation)

  # Get PC values for each coord
  pcavals <- data.frame(terra::extract(pcaRast, coords, ID = FALSE))
  colnames(pcavals) <- colnames(xpc)

  # Rescale var loadings with individual loadings so they fit in the plot nicely
  scldat <- min(
    (max(pcavals[, y], na.rm = TRUE) - min(pcavals[, y], na.rm = TRUE) / (max(varpc[, y], na.rm = TRUE) - min(varpc[, y], na.rm = TRUE))),
    (max(pcavals[, x], na.rm = TRUE) - min(pcavals[, x], na.rm = TRUE) / (max(varpc[, x], na.rm = TRUE) - min(varpc[, x], na.rm = TRUE)))
  )

  # Additionally use a constant scale val (scl) to shrink the final vectors (again for plotting nicely)
  varpc <- data.frame(varpc,
    v1 = scl * scldat * varpc[, x],
    v2 = scl * scldat * varpc[, y]
  )

  # GET RGB VALS FOR EACH COORD----------------------------------------------------------------------------------------

  pcavalsRGB <- data.frame(terra::extract(pcaRastRGB, coords, ID = FALSE))
  colnames(pcavalsRGB) <- colnames(xpc)

  # Create vector of RGB colors for plotting
  pcacols <- apply(pcavalsRGB, 1, create_rgb_vec)

  # GET RGB VALS FOR ENTIRE RASTER-------------------------------------------------------------------------------------

  # Get sample
  s <- sample(1:terra::ncell(pcaRast), 10000)

  # Get all PC values from raster and remove NAs
  rastvals <- data.frame(terra::values(pcaRast))[s, ]
  colnames(rastvals) <- colnames(xpc)
  rastvals <- rastvals[stats::complete.cases(rastvals), ]

  # Get all RGB values from raster and remove NAs
  rastvalsRGB <- data.frame(terra::values(pcaRastRGB))[s, ]
  colnames(rastvalsRGB) <- colnames(rastvals)
  rastvalsRGB <- rastvalsRGB[stats::complete.cases(rastvalsRGB), ]

  # Create vector of RGB colors for plotting
  rastpcacols <- apply(rastvalsRGB, 1, create_rgb_vec)

  # FINAL PLOT----------------------------------------------------------------------------------------------------------

  # Build base plot
  # Plot points colored by RGB with variable vectors
  plot <- ggplot2::ggplot() +

    # Create axes that cross through origin
    {
      if (display_axes) ggplot2::geom_hline(yintercept = 0, size = 0.2, col = "gray")
    } +
    {
      if (display_axes) ggplot2::geom_vline(xintercept = 0, size = 0.2, col = "gray")
    } +

    # Plot points from entire raster
    ggplot2::geom_point(data = rastvals, ggplot2::aes_string(x = x, y = y), col = rastpcacols, size = 4, alpha = 0.02) +

    # Plot coord values
    ggplot2::geom_point(data = pcavals, ggplot2::aes_string(x = x, y = y), fill = pcacols, col = "black", pch = 21, size = 3) +

    # Plot variable vectors
    ggplot2::geom_text(data = varpc, ggplot2::aes(x = v1, y = v2, label = varnames), size = 4, vjust = 1) +
    ggplot2::geom_segment(data = varpc, ggplot2::aes(x = 0, y = 0, xend = v1, yend = v2), arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))) +

    # Plot formatting
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      aspect.ratio = 1
    )

  # Build plot without PC axes displayed
  if (display_axes == FALSE) {
    plot <- plot +
      # Remove axes
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  }

  # Plot
  print(plot)
}


#' Helper function to create rgb vector
#'
#' @export
#' @noRd
create_rgb_vec <- function(vec) {
  if (any(is.na(vec))) x <- NA else x <- rgb(vec[1], vec[2], vec[3], maxColorValue = 255)
  return(x)
}


#' Scale a raster stack from 0 to 255
#'
#' @param s RasterStack
#'
#' @noRd
#' @export
stack_to_rgb <- function(s) {
  stack_list <- as.list(s)
  new_stack <- terra::rast(purrr::map(stack_list, raster_to_rgb))
  return(new_stack)
}


#' Scale raster from 0 to 255
#'
#' @param r SpatRast
#'
#' @noRd
#' @export
raster_to_rgb <- function(r) {
  rmax <- terra::minmax(r)["max", ]
  rmin <- terra::minmax(r)["min", ]
  if ((rmax - rmin) == 0) {
    r[] <- 255
  } else {
    r <- (r - rmin) / (rmax - rmin) * 255
  }
  return(r)
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
gdm_coeffs <- function(gdm_model) {
  # Vector to store coefficient sums
  coefSums <- c()

  # Sum coefficients for each predictor (each has 3 splines which you sum)
  for (i in 1:length(gdm_model$predictors)) {
    # Create starting index (* 3 b/c there are 3 splines and - 2 to get the first index)
    j <- (i * 3) - 2

    # Subset out three splines and sum them (j = first index, j + 2 = last index)
    coefSums[i] <- sum(gdm_model$coefficients[j:(j + 2)])
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
gdm_df <- function(gdm_result) {
  coeff_df <- gdm_coeffs(gdm_result$model)
  if (!is.null(gdm_result$pvalues)) coeff_df$p <- gdm_result$pvalues
  return(coeff_df)
}


#' Create `gt` table of GDM results
#'
#' @param gdm_result output of \link[algatr]{gdm_run} or \link[algatr]{gdm_do_everything} or a GDM model object
#' @param digits number of digits to include (defaults to 2)
#'
#' @family GDM functions
#'
#' @return An object of class `gt_tbl`
#' @export
gdm_table <- function(gdm_result, digits = 2, summary_stats = TRUE, footnote = TRUE) {

  # Format GDM model output as a list
  if (inherits(gdm_result, "gdm")) gdm_result <- list(model = gdm_result)
  
  gdm_df <- gdm_df(gdm_result)

  d <- max(abs(min(gdm_df$coefficient, na.rm = TRUE)), abs(max(gdm_df$coefficient)))

  suppressWarnings({
    tbl <- gdm_df %>%
      gt::gt() %>%
      gtExtras::gt_hulk_col_numeric(coefficient, trim = TRUE, domain = c(-d, d)) %>%
      gt::sub_missing(missing_text = "")

    if (summary_stats) {
      stat_names <- c("% Explained:")
      stats <- as.numeric(gdm_result$model$explained)
      gdm_df <- gdm_df %>%
        rbind(purrr::map2_dfr(.x = stat_names, .y = stats, .f = make_stat_vec, gdm_df)) %>%
        dplyr::mutate(dplyr::across(-c(predictor), as.numeric))


      tbl <- gdm_df %>%
        gt::gt() %>%
        gtExtras::gt_hulk_col_numeric(coefficient, trim = TRUE, domain = c(-d, d)) %>%
        gt::sub_missing(missing_text = "") %>%
        gt::tab_row_group(label = NA, id = "model", rows = which(!(gdm_df$predictor %in% stat_names))) %>%
        gtExtras::gt_highlight_rows(rows = which(gdm_df$predictor %in% stat_names), fill = "white") %>%
        gt::tab_style(
          style = list(
            gt::cell_borders(sides = "top", color = "white"),
            gt::cell_text(align = "left")
          ),
          locations = gt::cells_body(rows = which(gdm_df$var %in% stat_names))
        )
    }

    if (footnote & summary_stats) {
      tbl <- tbl %>% gt::tab_footnote(
        footnote = "The percentage of null deviance explained by the fitted GDM model.",
        locations = gt::cells_body(rows = which(gdm_df$predictor == "% Explained:"), columns = coefficient),
        placement = "right"
      )
    }


    if (!is.null(digits)) tbl <- tbl %>% gt::fmt_number(columns = -predictor, decimals = 2)
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
scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


#' Generate a Variable Importance Table for GDM Models
#'
#' This function generates a table displaying the variable importance for Generalized Dissimilarity Models (GDM). 
#' It can take either a `gdmData` object or a precomputed variable importance object and outputs a formatted table.
#'
#' @param varimp a `gdmData` object or a variable importance object created by running \link[gdm]{gdm.varImp}.
#' @param digits number of digits to include (defaults to 2)
#' @param summary_stats whether to add summary statistics to bottom of table (defaults to TRUE).
#' @param nPerm number of permutations to use if `varimp` is a `gdmData` object. Default is 50.
#' @param geo whether to include geographic distance in the GDM model. Default is TRUE.
#'
#' @return A `gt` table object displaying the variable importance.
#'
#' @export
gdm_varimp_table <- function(varimp, digits = 2, summary_stats = TRUE, nPerm = 50, geo = TRUE) {

  # If varimp is a gdmData object, run gdm.varImp
  if (inherits(varimp, "gdmData")) {
    varimp <- gdm::gdm.varImp(varimp, predSelect = FALSE, nPerm = nPerm, geo = geo)
  }

  # Make df
  df <- 
    varimp[-1] %>% 
    purrr::imap(\(x, i) {
      colnames(x) <- i
      x$Predictor <- row.names(x)
      return(x)
    }) %>% 
    purrr::reduce(dplyr::left_join, by = "Predictor") %>% 
    dplyr::select(Predictor, dplyr::everything())

  varstats <- varimp[[1]]

  # Round decimal places based on digits
  if (digits) df$`Predictor Importance` <- round(df$`Predictor Importance`, digits)
  d <- max(abs(min(df$`Predictor Importance`)), abs(max(df$`Predictor Importance`)))

  # Build table
  suppressWarnings({
    tbl <- 
      df %>% 
      gt::gt() %>% 
      gtExtras::gt_hulk_col_numeric(`Predictor Importance`, trim = TRUE, domain = c(-d, d)) %>% 
      gt::sub_missing(missing_text = "")

    # Add summary stats to bottom of table
    if (summary_stats) {
      stat_names <- c("Model deviance:", "Percent deviance explained:", "Model p-value:", "Fitted permutations:")
      stats <- varstats[,"All predictors"]
      df <- 
        df %>% 
        rbind(purrr::map2_dfr(.x = stat_names, .y = stats, .f = make_stat_vec, df)) %>% 
        dplyr::mutate(dplyr::across(-c(Predictor), as.numeric))

      tbl <- 
        df %>% 
        gt::gt() %>% 
        gtExtras::gt_hulk_col_numeric(`Predictor Importance`, trim = TRUE, domain = c(-d, d)) %>% 
        gt::sub_missing(missing_text = "") %>% 
        gt::tab_row_group(label = NA, id = "model", rows = which(!(df$Predictor %in% stat_names))) %>% 
        gtExtras::gt_highlight_rows(rows = which(df$Predictor %in% stat_names), fill = "white") %>% 
        gt::tab_style(
          style = list(
            gt::cell_borders(sides = "top", color = "white"),
            gt::cell_text(align = "left"),
            "padding-top:2px;padding-bottom:2px;"
          ),
          locations = gt::cells_body(rows = which(df$Predictor %in% stat_names))
        )
    }

    if (!is.null(digits)) tbl <- tbl %>% gt::fmt_number(columns = -Predictor, decimals = 2)
  })

  tbl
}
