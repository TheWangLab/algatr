#' MMRR.varSel performs MMRR with backward elimination variable selection
#' @param Y is a dependent distance matrix
#' @param X is a list of independent distance matrices (with optional names)
#' @param nperm is the number of permutations to be used in significance tests. Default = 999.
#' @param stdz if TRUE then matrices will be standardized. Default = TRUE.

MMRR.varSel <- function(Y, X, nperm = 999, stdz = TRUE){
  # Fit full model
  mmrr.model <- MMRR(Y, X, nperm = nperm, stdz = stdz)
  pvals <- mmrr.model$tpvalue[-1] # Remove intercept p-value
  
  # Eliminate variable with highest p-value, re-fit, and continue until only significant variables remain
  while((max(pvals) > 0.05) & (length(pvals) > 1)){
    print(pvals)
    rem.var <- which(pvals == max(pvals))
    X <- X[-rem.var]
    mmrr.model <- MMRR(Y, X, nperm = nperm, stdz = stdz)
    pvals <- mmrr.model$tpvalue[-1]
  }
  
  # Repetitive test, but just to be sure
  if((length(pvals) == 1) & all(pvals > 0.05)){warning("No significant variable combo found, returning NULL object"); mmrr.model <- NULL}
  
  return(mmrr.model)
}

#' MMRR performs Multiple Matrix Regression with Randomization analysis
#' @param Y is a dependent distance matrix
#' @param X is a list of independent distance matrices (with optional names)
#' @param nperm is the number of permutations to be used in significance tests. Default = 999.
#' @param stdz if TRUE then matrices will be standardized. Default = TRUE.

MMRR <- function(Y, X, nperm = 999, stdz = TRUE){
  # Compute regression coefficients and test statistics
  nrowsY <- nrow(Y)
  y <- unfold(Y, stdz)
  if(is.null(names(X)))names(X) <- paste("X", 1:length(X),sep="")
  Xmats <- sapply(X, unfold, stdz = stdz)
  fit <- lm(y ~ Xmats)
  coeffs <- fit$coefficients
  summ <- summary(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  tprob <- rep(1,length(tstat))
  Fprob <- 1
  # Get conf interval
  conf_df <- confint(fit, names(fit$coefficients), level = 0.90)
  rownames(conf_df) <-  c("Intercept", names(X))
  
  # Perform permutations
  for(i in 1:nperm){
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm, stdz)
    fit <- lm(yperm ~ Xmats)
    summ <- summary(fit)
    Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
    tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
  }
  
  # Return values
  tp <- tprob/(nperm+1)
  Fp <- Fprob/(nperm+1)
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
  if(stdz == TRUE) x <- scale(x, center = TRUE, scale = TRUE) 
  return(x)
}


#' Make dist matrix from vector
#'
#' @param vec vector
#'
#' @return
#' @export
#'
#' @examples
vec_to_dist <- function(vec){
  vec_dist <- as.matrix(dist(vec, diag = TRUE, upper = TRUE))
  return(vec_dist)
}

mmrr_doEverything <- function(gendist, coords, envlayers, model = "best", nperm = 999){
  
  # Make env
  env <- raster::extract(envlayers, coords)
  env <- as_tibble(env)
  
  # get geodist
  geodist <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
  
  # Make env dist matrix
  Xdist <- map(env, vec_to_dist)
  Xdist[["geodist"]] <- geodist
  
  # make geodist mat
  Ydist <- as.matrix(gendist)
  
  # Run MMRR with var selection
  if(model == "best"){mod <- MMRR.varSel(Ydist, Xdist, nperm = nperm)}
  if(model == "full"){mod <- MMRR(Ydist, Xdist, nperm = nperm)}
  
  # If NULL, exit with NULL
  if(is.null(mod)){warning("model is NULL, returning NULL object"); return(NULL)}
    
  # Make nice df
  coeff_df <- data.frame(coeff = mod$coefficients, p = mod$tpvalue)
  coeff_df$var <- rownames(coeff_df)
  ci_df <- data.frame(mod$conf_df)
  ci_df$var <- rownames(ci_df)
  coeff_df <- merge(coeff_df, ci_df, by = "var")
  rownames(coeff_df) <- NULL
  
  # Make results list
  results <- list(coeff_df = coeff_df,
                  mod = mod)
  
  # Return
  return(results)
  
}
