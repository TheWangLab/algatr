
#' Run RDA
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = loci)
#' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
#' @param coords dataframe with coordinates (only needed if correctGEO = TRUE) or if env is a Raster* from which values should be extracted
#' @param correctGEO whether to condition on geographic coordinates
#' @param correctPC whether to condition on PCs from PCA of genotypes
#' @param alpha significance level to use to identify loci
#' @param padj_method correction method supplied to \code{p.adjust} (can be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param nPC number of PCs to use if correctPC = TRUE
#'
#' @inheritParams vegan::ordiR2step
#' @inheritParams stats::p.adjust
#'
#' @return
#' @export
#'
#' @examples
rda_do_everything <- function(gen, env, coords = NULL, correctGEO = FALSE, correctPC = FALSE, alpha = 0.05, padj_method = "fdr", nPC = 3, Pin = 0.05, R2permutations = 1000, R2scope = T){

  # Modify environmental data --------------------------------------------------------------------------------------------------

  # extract environmental data if env is a RasterStack
  if(inherits(env, "Raster") | inherits(env, "RasterStack")) env <- raster::extract(env, coords)

  # Standardize environmental variables
  env <- scale(env, center = TRUE, scale = TRUE)
  env <- data.frame(env)

  # Rename coords
  colnames(coords) <- c("x", "y")

  # Check for NAs
  if(any(is.na(env))){
    warning("NA values found in env data, removing rows with NAs for RDA")
    gen <- gen[complete.cases(env),]
    coords <- coords[complete.cases(env),]
    # NOTE: this must be last
    env <- env[complete.cases(env),]
  }

  if(any(is.na(gen))){
    warning("NA values found in gen data, removing rows with NAs for RDA")
    env <- env[complete.cases(gen),]
    coords <- coords[complete.cases(gen),]
    # NOTE: this must be last
    gen <- gen[complete.cases(gen),]
  }


  # Running RDA ----------------------------------------------------------------------------------------------------------------

  # Run model
  mod <- rda_run(gen, env, coords, correctGEO, correctPC)

  # Identify candidate loci ----------------------------------------------------------------------------------------------------

  # Running with all axes
  rda_sig <- rda_getloci(mod, gen, naxes = NULL, padj_method = padj_method, alpha = alpha)

  # Return NULL if no loci are found
  if(is.null(rda_sig)) return(NULL)

  rda_loci <- rda_sig[["rda_loci"]]
  pvalues <- rda_sig[["pvalues"]]

  # Identify environmental associations
  rda_gen <- gen[,rda_loci]
  lm_df <- rda_lm_test(rda_gen, env)

  # Plot results ---------------------------------------------------------------------------------------------------------------

  # plot all axes
  rda_plot(mod, rda_loci, pvalues, biplot_axes = NULL, alpha = alpha)

  # Compile results ------------------------------------------------------------------------------------------------------------

  results <- list(rda_loci = rda_loci,
                  lm_df = lm_df,
                  rda_mod = mod,
                  pvalues = pvalues)

  return(results)
}


#' Run RDA
#'
#' @inheritParams rda_doEverything
#'
#' @return
#' @export
#'
#' @examples
rda_run <- function(gen, env, coords = NULL, correctGEO = FALSE, correctPC = FALSE, nPC = 3){

  if(!correctPC & !correctGEO){
    moddf <- data.frame(env)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse='+')))
  }

  if(correctPC & !correctGEO){
    pcres <- prcomp(gen)
    pc <-  pcres$x[,1:3]
    moddf <- data.frame(env, pc)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = '+'), "+ Condition(" , paste(colnames(pc), collapse = '+'), ")"))
  }

  if(!correctPC & correctGEO){
    if(is.null(coords)) stop("coordinates must be provided if correctGEO is TRUE")
    moddf <- data.frame(env, coords)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = '+'), "+ Condition(x + y)"))
  }

  if(correctPC & correctGEO){
    if(is.null(coords)) stop("coordinates must be provided if correctGEO is TRUE")
    pcres <- prcomp(gen)
    pc <-  pcres$x[,1:3]
    moddf <- data.frame(env, coords, pc)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = '+'), "+ Condition(" , paste(colnames(pc), collapse = '+'), "+ x + y)"))
  }


  mod <- vegan::rda(f, data = moddf)

  return(mod)

}



#' @export
rda_getloci <- function(mod, gen = NULL, naxes = NULL, padj_method = "fdr", alpha = 0.05){
  # Running the function with all axes
  # TODO: MODIFY THIS TO TAKE DIFFERENT NUMBERS OF AXES (AUTOMATICALLY CHOSEN USING SCREE or VAR explained)
  if(is.null(naxes)) naxes <- ncol(mod$CCA$v)
  rdadapt_env <- rdadapt(mod, naxes)

  # P-values threshold after FDR correction (different from Capblancq & Forester 2021)
  pvalues <- p.adjust(rdadapt_env$p.values, method = padj_method)

  # Capblancq include a step where they only take pvalues with highest loading for each contig to deal with LD (not applied here)
  # NOTE: I Think this filtering step should occur before (e.g. only one loci per LD block, but you know which comes from where)

  ## Identifying the loci that are below the p-value threshold
  #Identify rda cand loci (P)
  rda_loci <- which(pvalues < alpha)
  if (length(rda_loci) == 0) {
    warning("No significant loci found, returning NULL object")
    return(NULL)
  }

  # if a gen df is provided column names (loci names) are returned, otherwise indices are returned
  # NOTE: if indices are returned make sure to take extra steps to confirm these align properly
  if(!is.null(gen)){
    rda_loci <- colnames(gen)[rda_loci]
    names(pvalues) <- colnames(gen)
  } else {
    names(pvalues) <- 1:length(pvalues)
    message("no gen object provided - returning indices of significant loci")
  }

  results <- list(rda_loci = rda_loci,
                  pvalues = pvalues,
                  rdadapt = rdadapt_env)

  return(results)
}

# Function to conduct a RDA based genome scan from Capblancq & Forester 2021
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd
# NOTE: GO THROUGH THIS CODE
#' @export
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- robust::covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue::qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

#' Test for associations
#'
#' @param rda_gen genotype matrix for candidate loci
#' @param env env dataframe
#'
#' @return
#' @export
#'
#' @examples
rda_lm_test <- function(rda_gen, env){
  # suppress messages used because of "new names" message that is annoying
  env_p <- suppressMessages(purrr::map_dfc(rda_gen, rda_lm_env_helper, env))
  rda_df <- data.frame(colnames(env), env_p)
  colnames(rda_df) <- c("env", colnames(rda_gen))
  return(rda_df)
}

#' Helper function for rda_lm_test
#'
#' @param locus
#' @param env
#'
#' @keywords internal
#'
#' @return
#' @export
#'
#' @examples
rda_lm_env_helper <- function(locus, env){
  env_p <- purrr::map_df(env, rda_lm_helper, locus)
  colnames(env_p) <- names(locus)
  return(env_p)
}

#' Helper function for rda_lm_test
#'
#' @param envvar
#' @param locus
#'
#' @keywords internal
#'
#' @return
#' @export
#'
#' @examples
rda_lm_helper <- function(envvar, locus){
  mod_df <- data.frame(locus, envvar)
  colnames(mod_df) <- c("locus", "envvar")
  mod <- lm(locus ~ envvar, data = mod_df)
  mod_sum <- summary(mod)
  pvalue <- mod_sum$coefficients[,"Pr(>|t|)"]["envvar"]
  return(pvalue)
}

#' @export
rda_plot <- function(mod, rda_loci, pvalues, biplot_axes = NULL, alpha = 0.05, manhattan = TRUE, biplot = TRUE){
  # Get naxes
  naxes <- ncol(mod$CCA$v)

  # Make and get tidy dataframes for plotting
  tidy_list <- rda_tidy(mod, rda_loci, gen, naxes)
  TAB_loci <- tidy_list[["TAB_loci"]]
  TAB_var <- tidy_list[["TAB_var"]]

  # Make RDA plots
  if(biplot){
    if(!is.null(biplot_axes)){
      if(is.vector(biplot_axes)) print(rda_biplot(TAB_loci, TAB_var, biplot_axes = biplot_axes))
      if(is.list(biplot_axes)) lapply(biplot_axes, function(x) {print(rda_biplot(TAB_loci, TAB_var, biplot_axes = x))})
    } else {
      cb <- combn(naxes, 2)
      apply(cb, 2, function(x) {print(rda_biplot(TAB_loci, TAB_var, biplot_axes = x))})
    }
  }

  # Make manhattan plot
  if(manhattan) print(rda_manhattan(TAB_loci, rda_loci, pvalues, alpha = alpha))
}


#' Make dataframe for ggplot from RDA results
#'
#' @param rda_mod rda model
#' @param rda_loci rda candidate loci
#' @param naxes number of axes to include
#'
#' @return
#' @export
#'
#' @examples
rda_tidy <- function(mod, rda_loci, gen = NULL, naxes){
  locus_scores <- vegan::scores(mod, choices = 1:naxes, display="species", scaling="none") # vegan references "species", here these are the loci
  TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)

  TAB_loci$type <- "Neutral"
  TAB_loci$type[TAB_loci$names %in% rda_loci] <- "Outliers"
  TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "Outliers"))
  TAB_var <- as.data.frame(vegan::scores(mod, choices=1:naxes, display="bp")) # pull the biplot scores

  tidy_list <- list(TAB_loci = TAB_loci, TAB_var = TAB_var)
  return(tidy_list)
}


#' Plot RDA biplot
#'
#' @param TAB_loci
#' @param TAB_var
#' @param biplot_axes
#'
#' @return
#' @export
#'
#' @examples
rda_biplot <- function(TAB_loci, TAB_var, biplot_axes = c(1,2)){

  # Select axes for plotting
  xax <- paste0("RDA",biplot_axes[1])
  yax <- paste0("RDA",biplot_axes[2])
  TAB_loci_sub <- TAB_loci[,c(xax, yax, "type")]
  colnames(TAB_loci_sub) <- c("x", "y", "type")
  TAB_var_sub <- TAB_var[,c(xax, yax)]
  colnames(TAB_var_sub) <- c("x", "y")

  ## Biplot of RDA loci and variables scores
  ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_point(data = TAB_loci_sub, aes(x=x*20, y=y*20, colour = type), size = 1.4) +
    scale_color_manual(values = c(rgb(0.7,0.7,0.7,0.1), "#F9A242FF")) +
    geom_segment(data = TAB_var_sub, aes(xend=x, yend=y, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_text(data = TAB_var_sub, aes(x=1.1*x, y=1.1*y, label = row.names(TAB_var_sub)), size = 2.5) +
    xlab(xax) + ylab(yax) +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11) +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
}

#' RDA manhattan plot
#'
#' @param TAB_loci
#' @param pvalues
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
rda_manhattan <- function(TAB_loci, rda_loci, pvalues, alpha = 0.05){

  ## Manhattan plot
  TAB_manhattan <- data.frame(pos = 1:nrow(TAB_loci),
                              pvalues = pvalues,
                              type = factor(TAB_loci$type, levels = c("Neutral", "Outliers")))
  TAB_manhattan <- TAB_manhattan[order(TAB_manhattan$pos),]

  ggplot(data = TAB_manhattan) +
    geom_point(aes(x=pos, y=-log10(pvalues), col = type), size=1.4) +
    scale_color_manual(values = c(rgb(0.7,0.7,0.7,0.5), "#F9A242FF", "#6B4596FF")) +
    xlab("Loci") + ylab("-log10(p)") +
    geom_hline(yintercept=-log10(alpha), linetype="dashed", color = "black", size=0.6) +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11) +
    theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
}

#' Make pretty table of p-values
#'
#' @param lm_df
#'
#' @return
#' @export
#'
#' @examples
lm_table <- function(lm_df){
  suppressWarnings(
    tbl <-lm_df %>%
      mutate(across(-env, round, 3)) %>%
      gt() %>%
      data_color(columns = -env, colors = scales::col_numeric(
        palette = c("#F9A242FF"),
        domain = c(0, 0.05)
      )
      )
  )
  tbl
}

