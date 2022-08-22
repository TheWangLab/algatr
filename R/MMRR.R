
#' MMRR function to do everything ()
#' TODO [EAC]: gendist matrix must range between 0 and 1 as with GDM?
#'
#' @param gendist matrix of genetic distances
#' @param coords dataframe with x and y coordinates (MUST BE CALLED X AND Y)
#' @param envlayers envlayers for mapping (MUST MATCH NAMES IN ENV DATAFRAME)
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "best"
#' @param n_perm number of permutations to use to calculate variable importance, only matters if model = "best" (defaults to 999)
#' @param env dataframe with environmental values for each coordinate, if not provided it will be calculated based on coords/envlayers
#'
#' @return
#' @export
#'
#' @examples
mmrr_do_everything <- function(gendist, coords, env = NULL, envlayers, model = "best", n_perm = 999){

  # If not provided, make env data frame from layers and coords
  if(is.null(env)){env <- raster::extract(envlayers, coords)}

  # Convert to tibble for map function
  env <- dplyr::as_tibble(env)

  # Make env dist matrix
  Xdist <- purrr::map(env, env_dist)
  Xdist[["geodist"]] <- geo_dist(coords)

  # Make geodist mat
  Ydist <- as.matrix(gendist)

  # Run MMRR with variable selection
  if(model == "best"){mod <- mmrr_var_sel(Ydist, Xdist, n_perm = n_perm)}
  if(model == "full"){mod <- MMRR(Ydist, Xdist, n_perm = n_perm)}

  # If NULL, exit with NULL
  if(is.null(mod)){warning("model is NULL, returning NULL object"); return(NULL)}

  # Make nice data frame
  coeff_df <- data.frame(coeff = mod$coefficients, p = mod$tpvalue)
  coeff_df$var <- rownames(coeff_df)
  ci_df <- data.frame(mod$conf_df)
  ci_df$var <- rownames(ci_df)
  coeff_df <- merge(coeff_df, ci_df, by = "var")
  rownames(coeff_df) <- NULL

  # Make results list
  results <- list(coeff_df = coeff_df,
                  mod = mod)

  return(results)

}

#' mmrr_var_sel performs MMRR with backward elimination variable selection
#' @param Y is a dependent distance matrix
#' @param X is a list of independent distance matrices (with optional names)
#' @param n_perm is the number of permutations to be used in significance tests. Default = 999.
#' @param stdz if TRUE then matrices will be standardized. Default = TRUE.

mmrr_var_sel <- function(Y, X, n_perm = 999, stdz = TRUE){
  # Fit full model
  mmrr.model <- MMRR(Y, X, n_perm = n_perm, stdz = stdz)
  pvals <- mmrr.model$tpvalue[-1] # Remove intercept p-value

  # Eliminate variable with highest p-value, re-fit, and continue until only significant variables remain
  while((max(pvals) > 0.05) & (length(pvals) > 1)){
    print(pvals)
    rem.var <- which(pvals == max(pvals))
    X <- X[-rem.var]
    mmrr.model <- MMRR(Y, X, n_perm = n_perm, stdz = stdz)
    pvals <- mmrr.model$tpvalue[-1]
  }

  # Repetitive test, but just to be sure
  if((length(pvals) == 1) & all(pvals > 0.05)){warning("No significant variable combo found, returning NULL object"); mmrr.model <- NULL}

  return(mmrr.model)
}

#' MMRR performs Multiple Matrix Regression with Randomization analysis
#' @param Y is a dependent distance matrix
#' @param X is a list of independent distance matrices (with optional names)
#' @param n_perm is the number of permutations to be used in significance tests. Default = 999.
#' @param stdz if TRUE then matrices will be standardized. Default = TRUE.

MMRR <- function(Y, X, n_perm = 999, stdz = TRUE){
  # Compute regression coefficients and test statistics
  nrowsY <- nrow(Y)
  y <- unfold(Y, stdz)
  if(is.null(names(X)))names(X) <- paste("X", 1:length(X),sep="")
  Xmats <- sapply(X, unfold, stdz = stdz)
  fit <- stats::lm(y ~ Xmats)
  coeffs <- fit$coefficients
  summ <- summary(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  tprob <- rep(1,length(tstat))
  Fprob <- 1

  # Get confidence interval
  conf_df <- stats::confint(fit, names(fit$coefficients), level = 0.90)
  rownames(conf_df) <-  c("Intercept", names(X))

  # Perform permutations
  for(i in 1:n_perm){
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm, stdz)
    fit <- stats::lm(yperm ~ Xmats)
    summ <- summary(fit)
    Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
    tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
  }

  # Return values
  tp <- tprob/(n_perm+1)
  Fp <- Fprob/(n_perm+1)
  names(r.squared) <- "r.squared"
  names(coeffs) <- c("Intercept", names(X))
  names(tstat) <- paste(c("Intercept", names(X)), "(t)", sep="")
  names(tp) <- paste(c("Intercept", names(X)), "(p)", sep="")
  names(Fstat) <- "F-statistic"
  names(Fp) <- "F p-value"

  return(list(r.squared = r.squared,
              coefficients = coeffs,
              tstatistic = tstat,
              tpvalue = tp,
              Fstatistic = Fstat,
              Fpvalue = Fp,
              conf_df = conf_df))
}

#' unfold converts the lower diagonal elements of a matrix into a vector
#' @param X is a distance matrix.
#' @param stdz if TRUE then matrices will be standardized. Default = TRUE.

unfold <- function(X, stdz = TRUE){
  x <- vector()
  for(i in 2:nrow(X)) x <- c(x, X[i, 1:i-1])
  if(stdz == TRUE) x <- raster::scale(x, center = TRUE, scale = TRUE)
  return(x)
}

