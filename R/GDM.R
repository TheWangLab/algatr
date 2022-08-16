
#' GDM function to do everything (fit model, get coefficients, make and save raster)
#'
#' @param gendist matrix of genetic distances (MUST RANGE BETWEEN 0 AND 1)
#' @param coords dataframe with x and y coordinates (MUST BE CALLED X AND Y) TODO: FIX THIS
#' @param env dataframe with environmental values for each coordinate, if not provided it will be calculated based on coords/envlayers
#' @param envlayers envlayers for mapping (MUST MATCH NAMES IN ENV DATAFRAME)
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "best"
#' @param alpha alpha level for variable selection (defaults to 0.05), only used if model = "best" TODO: ADD BETTER DESCRIPTOR
#' @param nPerm number of permutations to use to calculate variable importance, only matters if model = "best" (defaults to 50)
#' @param scale whether to scale genetic distance data from 0 to 1 (defaults to FALSE) TODO: ADD STOP
#' @param plot_vars whether to create variable vector loading plot (defaults to TRUE)
#'
#' @return list with final model, predictor coefficients, and PCA RGB map
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_do_everything <- function(gendist, coords, env = NULL, envlayers = NULL, model = "best", alpha = 0.05, nPerm = 50, scale = FALSE, plot_vars = TRUE){

  # TODO: insert a CHECK here to match sample IDs between gendist and coords, otherwise insert warning
  # to remind user to make sure before running anything

  # If not provided, make env data frame from layers and coords
  if(is.null(env)){env <- raster::extract(envlayers, coords)}

  # Run model with all defaults
  mod <- gdm_run(gendist, coords, env, model = model, alpha = alpha, nPerm = nPerm, scale = scale)

  # If mod is null, exit
  if(is.null(mod)){warning("GDM model is NULL, returning NULL object"); return(NULL)}

  # Get coefficients from models
  predictors <- gdm_coeffs(mod)

  # Create and plot map
  map <- gdm_map(mod, envlayers, coords, plot_vars = plot_vars)

  # Plot isplines
  gdm_plot_isplines(mod)

  # Create list to store results
  results <- list()
  # Add model
  results[["model"]] <- mod
  # Add predictors
  results[["predictors"]] <- predictors
  # Add raster(s)
  results[["rast"]] <- map

  return(results)
}


#' Run GDM and return model object
#'
#' @param gendist matrix of genetic distances (MUST RANGE BETWEEN 0 AND 1)
#' @param coords dataframe with x and y coordinates
#' @param env dataframe with environmental values for each coordinate
#' @param model whether to compute the full model ("full") or the best model based on variable selection steps ("best")
#' @param alpha alpha level for variable selection (defaults to 0.05), only matters if model = "best"
#' @param nPerm number of permutations to use to calculate variable importance, only matters if model = "best"
#' @param scale scale genetic distance data from 0 to 1
#'
#' @family GDM functions
#'
#' @return GDM model
#' @export
#'
#' @examples
gdm_run <- function(gendist, coords, env, model = "best", alpha = 0.05, nPerm = 50, scale = FALSE){

  # FORMAT DATA ---------------------------------------------------------------------------------------------------
  # TODO: CHECK SAMPLE IDS between gendist and coords; remove sample IDs from coords and continue

  # Rename coords
  coords <- dplyr::as_tibble(coords)
  colnames(coords) <- c("x","y")

  # Scale genetic distance data from 0 to 1
  if(scale){gendist <- scale01(gendist)}

  # Vector of sites (for individual-based sampling, this is just assigning 1 site to each individual)
  # TODO: should there be a warning here specifying that ind-level sampling is assumed?
  site <- 1:nrow(gendist)

  # Bind vector of sites with gen distances
  gdmGen <- cbind(site, gendist)

  # Create data frame of predictor variables
  gdmPred <- data.frame(site = site,
                        x = coords$x,
                        y = coords$y,
                        env)

  # Format data for GDM
  gdmData <- gdm::formatsitepair(gdmGen, 3, XColumn = "x", YColumn = "y", siteColumn = "site", predData = gdmPred)

  # RUN GDM -------------------------------------------------------------------------------------------------------

  # If model = "full", the final GDM model is just the full model
  if(model == "full"){
    # Remove any remaining incomplete cases
    cc <- stats::complete.cases(gdmData)
    if(!all(cc)){gdmData <- gdmData[cc, ]; warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))}

    # Run GDM with all predictors
    gdm_model_final <- gdm::gdm(gdmData, geo = TRUE)
  }

  # If model = "best", go through variable selection procedure
  if(model == "best"){
    # Add distance matrix separately
    geodist <- geo_dist(coords[,c("x","y")])
    gdmDist <- cbind(site, geodist)
    gdmData <- gdm::formatsitepair(gdmData, 4, predData = gdmPred, siteColumn = "site", distPreds = list(geodist = as.matrix(gdmDist)))

    # Remove any remaining incomplete cases
    cc <- stats::complete.cases(gdmData)
    if(!all(cc)){gdmData <- gdmData[cc, ]; warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))}

    # Get subset of variables for final model
    finalvars <- gdm_var_select(gdmData, alpha = alpha, nPerm = nPerm)

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

    return(gdm_model_final)

  }

  return(gdm_model_final)
}



#' Get best set of variables from a GDM model
#'
#' @param gdmData data formatted using GDM package
#' @param alpha alpha level for determining variable significance
#' @param nPerm number of permutations to run for variable testing
#'
#' @return
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_var_select <- function(gdmData, alpha = 0.05, nPerm = 10){

  # Check var importance/significance (THIS STEP CAN TAKE A WHILE)
  # TODO: why geo default is false?
  vars <- gdm::gdm.varImp(gdmData,
                     geo = FALSE,
                     splines = NULL,
                     nPerm = nPerm)

  # Get pvalues from variable selection model
  pvalues <- vars[[3]]

  # Identify which cells have p-values lower than alpha and not NA
  # NB: NA occurs when variables are removed during model testing
  cond <- (pvalues < alpha) | is.na(pvalues)

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

  # If the geodist matrix (matrix_1) is a significant variable, add geo to the list of vars
  if("matrix_1" %in% finalvars) {
    # Remove matrix
    finalvars <- finalvars[which(finalvars != "matrix_1")]
    # Add geo
    finalvars <- c("geo", finalvars)
  }

  return(finalvars)
}



#' Make map from model
#'
#' @param gdm_model gdm model
#' @param envlayers stack of raster layers (NAMES MUST CORRESPOND WITH GDM MODEL)
#' @param plot_vars whether to create PCA plot to help in variable and map interpretation
#' @param coords
#' @param scl # TODO: what is this?
#' @param plot
#'
#' @return GDM RGB map
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_map <- function(gdm_model, envlayers, coords, plot_vars = TRUE, scl = 1, plot = TRUE){

  # CHECK that all of the model variables are included in the stack of environmental layers
  # Create list of environmental predictors (everything but Geographic)
  check_geo <- gdm_model$predictors == "Geographic"
  if(any(check_geo)){model_vars <- gdm_model$predictors[-which(check_geo)]} else {model_vars <- gdm_model$predictors}

  # Check that model variables are included in names of envlayers
  var_check <- model_vars %in% names(envlayers)

  # Print error with missing layers
  if(!all(var_check)){stop(paste("missing model variable(s) from raster stack:",  model_vars[!var_check]))}

  # Subset envlayers to only include variables in final model
  envlayers_sub <- raster::subset(envlayers, model_vars)


  # CREATE MAP ----------------------------------------------------------------------------------------------------

  # Transform GIS layers
  rastTrans <- gdm::gdm.transform(gdm_model, envlayers_sub)

  # Remove NA values
  rastDat <- na.omit(raster::getValues(rastTrans))

  # Run PCA
  pcaSamp <- stats::prcomp(rastDat)

  # Count number of layers
  n_layers <- raster::nlayers(rastTrans)
  # Max number of layers to plot is 3, so adjust n_layers accordingly
  if(n_layers > 3){n_layers <- 3}

  # Make PCA raster
  pcaRast <- raster::predict(rastTrans, pcaSamp, index=1:n_layers)

  # Make copy of PCA raster to transform into RGB values
  pcaRastRGB <- pcaRast

  # Scale rasters to get colors (each layer will correspond with R, G, or B in the final plot)
  for(layer in 1:n_layers){
    # If all values are 0 assign value as white (255)
    if(pcaRastRGB[[layer]]@data@max == 0){pcaRastRGB[[layer]] <- 255; next}
    pcaRastRGB[[layer]] <- (pcaRastRGB[[layer]]-pcaRastRGB[[layer]]@data@min) / (pcaRastRGB[[layer]]@data@max-pcaRastRGB[[layer]]@data@min)*255
  }

  # If there are fewer than 3 n_layers (e.g., <3 variables), the RGB plot won't work (because there isn't an R, G, and B)
  # To get around this, create a blank raster (i.e., a white raster), and add it to the stack
  # TODO: !COME BACK TO THIS!
  if(n_layers < 3){
    warning("Fewer than three non-zero coefficients provided, adding white substitute layers to RGB plot")
    # Create white raster by multiplying a layer of pcaRast by 0 and adding 255
    white_raster <- pcaRastRGB[[1]]*0 + 255
  }

  # If n_layers = 2, you end up making a bivariate map
  if(n_layers == 2){pcaRastRGB <- raster::stack(pcaRastRGB, white_raster)}
  # If n_layers = 1, you end up making a univariate map
  if(n_layers == 1){pcaRastRGB <- raster::stack(pcaRastRGB, white_raster, white_raster)}

  # Plot raster
  # TODO: points() in terra?
  if(plot){plotRGB(pcaRastRGB, r = 1, g = 2, b = 3)}
  if(!is.null(coords)) points(coords, cex = 1.5)

  # Plot variable vectors
  if(plot_vars & (n_layers == 3)){
    # TODO: FIX SO THAT THIS WORKS IF NO COORDS ARE PROVIDED
    gdm_plot_vars(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = scl)
  }
  if(plot_vars & (n_layers != 3)){
    warning("variable vector plot is not available for model with fewer than 3 final variables, skipping...")
  }

  s <- list(rastTrans, pcaRastRGB)
  names(s) <- c("rastTrans", "pcaRastRGB")
  return(s)

}


#' Plot isplines for each variable
#'
#' @param gdm_model GDM model
#'
#' @return plot for each ispline
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_plot_isplines <- function(gdm_model){
  gdm_model_splineDat <- isplineExtract(gdm_model)

  par(mfrow=c(1,ncol(gdm_model_splineDat$x)), pty = "s")

  for(i in 1:ncol(gdm_model_splineDat$x)){
    plot(gdm_model_splineDat$x[,i],
         gdm_model_splineDat$y[,i],
         cex.lab = 1.5,
         lwd = 2,
         type = "l",
         col = "black",
         ylim = c(0, max(gdm_model_splineDat$y)),
         xlab = colnames(gdm_model_splineDat$x)[i],
         ylab = "Partial Regression Distance")

  }
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
#'
#' @return GDM PCA plot
#'
#' @family GDM functions
#'
#' @export
#'
#' @examples
gdm_plot_vars <- function(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = 0.5){

  # Confirm there are exactly 3 axes
  if(nlayers(pcaRastRGB) > 3){stop("Only three PC layers (RGB) can be used for creating the variable plot (too many provided)")}
  if(nlayers(pcaRastRGB) < 3){stop("Need exactly three PC layers (RGB) for creating the variable plot (too few provided)")}

  # GET PCA DATA ----------------------------------------------------------------------------------------------------

  # Make data frame from PC results
  xpc <- data.frame(pcaSamp$x[,1:3])

  # Get variable rotations
  varpc <- data.frame(varnames = rownames(pcaSamp$rotation), pcaSamp$rotation)

  # Get PC values for each coord
  pcavals <- data.frame(raster::extract(pcaRast, coords))
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

  pcavalsRBG <- data.frame(raster::extract(pcaRastRGB, coords))
  colnames(pcavalsRBG) <- colnames(xpc)

  # Create vector of RGB colors for plotting
  pcacols <- c()
  for(i in 1:nrow(pcavalsRBG)){
    if(any(is.na(pcavalsRBG[i,]))){
      pcacols[i] <- NA
    }
    else {
      # r=1, g=2, b=3 (TODO: FIX TO MAKE THIS FLEXIBLE)
      pcacols[i] <- rgb(pcavalsRBG[,1][i], pcavalsRBG[,2][i], pcavalsRBG[,3][i], maxColorValue = 255)
    }

  }


  # GET RGB VALS FOR ENTIRE RASTER-------------------------------------------------------------------------------------

  # Get sample
  s <- sample(1:ncell(pcaRast), 10000)

  # Get all PC values from raster and remove NAs
  rastvals <- data.frame(values(pcaRast))[s,]
  colnames(rastvals) <- colnames(xpc)
  rastvals <- rastvals[complete.cases(rastvals),]

  # Get all RGB values from raster and remove NAs
  rastvalsRGB <- data.frame(values(pcaRastRGB))[s,]
  colnames(rastvalsRGB) <- colnames(rastvals)
  rastvalsRGB <- rastvalsRGB[complete.cases(rastvalsRGB),]

  # Create vector of RGB colors for plotting
  rastpcacols <- c()
  for(i in 1:nrow(rastvalsRGB)){
    if(any(is.na(rastvalsRGB[i,]))){
      rastpcacols[i] <- NA
    } else {
      # r=1, g=2, b=3 (TODO: FIX TO MAKE THIS FLEXIBLE)
      rastpcacols[i] <- rgb(rastvalsRGB[,1][i], rastvalsRGB[,2][i], rastvalsRGB[,3][i], maxColorValue = 255)
    }

  }

  # FINAL PLOT----------------------------------------------------------------------------------------------------------

  # Plot points colored by RGB with variable vectors
  plot <- ggplot() +

    # Create axes that cross through origin
    geom_hline(yintercept = 0, size=0.2, col = "gray") +
    geom_vline(xintercept = 0, size=0.2, col = "gray") +

    # Plot points from entire raster
    geom_point(data = rastvals, aes_string(x=x, y=y), col = rastpcacols, size = 4, alpha = 0.02) +

    # Plot coord values
    geom_point(data = pcavals, aes_string(x=x, y=y), fill = pcacols, col = "black", pch = 21, size = 3) +

    # Plot variable vectors
    geom_text(data = varpc, aes(x=v1, y=v2, label=varnames), size = 4, vjust=1) +
    geom_segment(data = varpc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm"))) +

    # Plot formatting
    coord_equal() +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(), aspect.ratio=1)

  # Plot
  print(plot)
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
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}


#' Get coefficients for each predictor
#'
#' @param gdm_model - model of type gdm
#'
#' @return data frame of coefficients for GDM
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

  # Add those values to a simple data frame
  coeffs <- data.frame(predictor = gdm_model$predictors, coefficient = coefSums)

  return(coeffs)
}


# TODO: ADD THESE IN LATER
stack_to_rgb <- function(s){
  new_stack <- s
  for(i in 1:nlayers(s)){new_stack[[i]] <- raster_to_rgb(s[[i]])}
  return(new_stack)
}

raster_to_rgb <- function(r){
  if(r@data@max == 0){r <- 255} else {r <- (r-r@data@min) / (r@data@max-r@data@min)*255}
}
