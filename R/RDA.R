
#' Run RDA
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = snps) or `vcfR` object
#' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
#' @param coords dataframe with coordinates (only needed if correctGEO = TRUE) or if env is a Raster* from which values should be extracted
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "best"
#' @param correctGEO whether to condition on geographic coordinates
#' @param correctPC whether to condition on PCs from PCA of genotypes
#' @param outlier_method method to determine outliers. Can either be "p" to use the p-value method from https://github.com/Capblancq/RDA-landscape-genomics or "z" to use the z-score based method from https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#' @param sig if `outlier_method = "p"`, the significance level to use to identify snps
#' @param padj_method if `outlier_method = "p"`, the correction method supplied to \code{p.adjust} (can be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param z if `outlier_method = "z"`, the number of standard deviations to use to identify snps
#' @param cortest whether to create table of correlations for snps and environmental variable
#' @param nPC number of PCs to use if correctPC = TRUE (defaults to 3), if set to "manual" a selection option with a terminal prompt will be provided
#' @param naxes number of RDA axes to use (defaults to "all" to use all axes), if set to "manual" a selection option with a terminal prompt will be given, otherwise can be any integer that is less than or equal to the total number of axes
#' @param Pin if `model = "best"`, limits of permutation P-values for adding (`Pin`) a term to the model, or dropping (`Pout`) from the model. Term is added if` P <= Pin`, and removed if `P > Pout` (see \link[vegan]{ordi2step})
#' @param R2permurations if `model = "best"`, number of permutations used in the estimation of adjusted R2 for cca using RsquareAdj (see \link[vegan]{ordi2step})
#' @param R2scope if `model = "best"`, use adjusted R2 as the stopping criterion: only models with lower adjusted R2 than scope are accepted (see \link[vegan]{ordi2step})
#' @param stdz whether to center and scale environmental data (defaults to TRUE)
#' @param impute function to use for imputation (defaults to `median`). NOTE: use extreme caution when using this form of simplistic imputation. We mainly provide this code for creating test datasets and highly discourage its use in analyses.
#'
#' @inheritParams vegan::ordiR2step
#'
#' @importFrom vegan rda
#'
#' @return list containing (1) outlier SNPs, (2) data frame with correlation test results; if `cortest = TRUE`: (3) the RDA model, (4) results from outlier analysis (output from \link[algatr]{rda_getoutliers}), (5) RDA R-Squared, (6) RDA ANOVA, (7) p-values if `outlier_method = "p"`,
#' @export
#'
#' @examples
rda_do_everything <- function(gen, env, coords = NULL, model = "best", correctGEO = FALSE, correctPC = FALSE,
                              outlier_method = "p", sig = 0.05, z = 3,
                              padj_method = "fdr", cortest = TRUE, nPC = 3, naxes = "all",
                              Pin = 0.05, R2permutations = 1000, R2scope = T, stdz = TRUE, impute = "median"){

  # Modify environmental data --------------------------------------------------------------------------------------------------

  # Extract environmental data if env is a RasterStack
  if(inherits(env, "Raster")) env <- raster::extract(env, coords)

  # Standardize environmental variables
  if(stdz) env <- scale(env, center = TRUE, scale = TRUE)
  env <- data.frame(env)

  # Rename coords
  colnames(coords) <- c("x", "y")

  # Modify genetic data -----------------------------------------------------

  # Convert vcf to dosage
  if(inherits(gen, "vcfR")) gen <- vcf_to_dosage(gen)

  # Perform imputation with warning
  if(any(is.na(gen))){
    gen <- simple_impute(gen, median)
    warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
  }

  # Check for NAs
  if(any(is.na(gen))){
    stop("NA values found in gen data")
  }

  if(any(is.na(env))){
    warning("NA values found in env data, removing rows with NAs for RDA")
    gen <- gen[complete.cases(env),]
    coords <- coords[complete.cases(env),]
    # NOTE: this must be last
    env <- env[complete.cases(env),]
  }

  # Running RDA ----------------------------------------------------------------------------------------------------------------

  # Run model
  mod <- rda_run(gen, env, coords,
                 model = model,
                 correctGEO = correctGEO,
                 correctPC = correctPC,
                 nPC = nPC,
                 Pin = Pin,
                 R2permutations = R2permutations,
                 R2scope = R2scope)

  # If NULL, exit
  if(is.null(mod)) stop("model is NULL")

  # get RSquared and run ANOVA
  mod_rsq <- vegan::RsquareAdj(mod)
  mod_aov <- anova(mod)

  # Identify candidate snps ----------------------------------------------------------------------------------------------------

  # Running with all axes
  rda_sig <- rda_getoutliers(mod, naxes = naxes, outlier_method = outlier_method, padj_method = padj_method, sig = sig)

  # Get SNPs
  rda_snps <- rda_sig$rda_snps

  # Summarize results ---------------------------------------------------------------------------------------------------------------

  # plot all axes
  if(any("pvalues" %in% names(rda_sig))) pvalues <- rda_sig[["pvalues"]] else pvalues <- NULL
  rda_plot(mod, rda_snps = rda_snps, pvalues = pvalues,  axes = "all", biplot_axes = NULL, sig = sig)

  # get correlations -----------------------------------------------------------------------------------------------------------
  rda_gen <- gen[,rda_snps]
  if(cortest) {
    cor_df <- rda_cor(rda_gen, env)
    print(rda_table(cor_df, top = TRUE, order = TRUE, nrow = 10))
  } else cor_df <- NULL

  # Compile results ------------------------------------------------------------------------------------------------------------

  results <- list(rda_snps = rda_snps,
                  cor_df = cor_df,
                  rda_mod = mod,
                  rda_outlier_test = rda_sig,
                  rsq = mod_rsq,
                  anova = mod_aov,
                  pvalues = pvalues)

  return(results)
}


#' Run RDA
#'
#' @inheritParams rda_doEverything
#'
#' @return RDA model
#' @export
#'
rda_run <- function(gen, env, coords = NULL, model = "full",
                    correctGEO = FALSE, correctPC = FALSE, nPC = 3,
                    Pin = 0.05, R2permutations = 1000, R2scope = T){

  if(!correctPC & !correctGEO){
    moddf <- data.frame(env)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse='+')))
  }

  if(correctPC & !correctGEO){
    pcres <- prcomp(gen)
    stats::screeplot(pcres, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
    if(nPC == "manual") nPC <- readline("Number of PC axes to retain:")
    pc <-  pcres$x[,1:nPC]
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
    stats::screeplot(pcres, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
    if(nPC == "manual") nPC <- readline("Number of PC axes to retain:")
    pc <- pcres$x[,1:3]
    moddf <- data.frame(env, coords, pc)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = '+'), "+ Condition(" , paste(colnames(pc), collapse = '+'), "+ x + y)"))
  }

  if(model == "best"){
    mod_full <- vegan::rda(f, data = moddf)
    mod_null <- vegan::rda(gen ~ 1,  data = moddf)
    mod <- vegan::ordiR2step(mod_null, mod_full, Pin = Pin, R2permutations = R2permutations, R2scope = R2scope)
    if(mod$call == mod_null$call) {mod <- NULL; warning("Best model is NULL model, returning NULL")}
  } else {
    mod <- vegan::rda(f, data = moddf)
  }

  return(mod)

}



#' Get significant outliers from RDA model
#'
#' @param plot whether to produce scree plot of RDA axes (defaults to TRUE)
#' @inheritParams rda_do_everything
#'
#' @return results from outlier tests. If `outlier_method = "p"`, a list of outlier snps, p-values, and results from rdadapt test (see Capblancq & Forester 2021; https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd). If `outlier_method = "z"`, a dataframe with outlier snp z-scores for each axes
#'
#' @export
rda_getoutliers <- function(mod, naxes = "all", outlier_method = "p", padj_method = "fdr", sig = 0.05, z = 3, plot = TRUE){
  # Running the function with all axes
  if(plot) stats::screeplot(mod, main = "Eigenvalues of constrained axes")
  if(naxes == "manual") naxes <- readline("Number of RDA axes to retain:")
  if(naxes == "all") naxes <- ncol(mod$CCA$v)

  if(outlier_method == "p" & naxes == 1) warning("Cannot compute p-values (outlier_method = \"p\") when the number of RDA axes is less than two, using the standard deviation based method (outlier_method = \"z\") instead")
  if(outlier_method == "p" & naxes != 1) results <- p_outlier_method(mod, naxes, sig, padj_method)
  if(outlier_method == "z" | naxes == 1) results <- z_outlier_method(mod, naxes, z)

  return(results)
}



#' Determine RDA outliers based on p-values
#' @inheritParams rda_getoutliers
#' @export
#' @noRd
#'
p_outlier_method <- function(mod, naxes, sig = 0.05, padj_method = "fdr"){
  rdadapt_env <- rdadapt(mod, naxes)

  # P-values threshold after FDR correction (different from Capblancq & Forester 2021)
  pvalues <- p.adjust(rdadapt_env$p.values, method = padj_method)

  # get snp names
  snp_names <- rownames(vegan::scores(mod, choices = naxes, display = "species"))

  # restore SNP names
  names(pvalues) <- snp_names

  # Capblancq include a step where they only take pvalues with highest loading for each contig to deal with LD (not applied here)
  # NOTE: I Think this filtering step should occur before (e.g. only one snps per LD block, but you know which comes from where)

  ## Identifying the snps that are below the p-value threshold
  #Identify rda cand snps (P)
  rda_snps <- snp_names[which(pvalues < sig)]
  if (length(rda_snps) == 0) {
    warning("No significant snps found, returning NULL object")
    return(NULL)
  }

  results <- list(rda_snps = rda_snps,
                  pvalues = pvalues,
                  rdadapt = rdadapt_env)

  return(results)
}

#' Determine RDA outliers based on Z-scores
#' @inheritParams rda_getoutliers
#' @export
#' @noRd
#'
z_outlier_method <- function(mod, naxes, z = 3){
  load.rda <- vegan::scores(mod, choices = naxes, display="species")

  results <- purrr::map_dfr(data.frame(1:ncol(load.rda)), z_outlier_helper, load.rda, z)

  return(results)
}

#' z_outlier_method helper function
#'
#' @export
#' @noRd
#'
z_outlier_helper <- function(axis, load.rda, z){
  x <- load.rda[,axis]
  out <- outliers(x, z)
  cand <- cbind.data.frame(names(out), rep(axis, times=length(out)), unname(out))
  colnames(cand) <- c("rda_snps", "axis", "loading")
  cand$rda_snps <- as.character(cand$rda_snps)
  return(cand)
}

#' Z outlier finder
#' from https://popgen.nescent.org/2018-03-27_RDA_GEA.html
#'
#' @export
#' @noRd
#'
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]               # snp names in these tails
}

# Function to conduct a RDA based genome scan from Capblancq & Forester 2021
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd
# TODO[to: EAC, from: APB]: GO THROUGH THIS CODE
#' @export
#' @noRd
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- robust::covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue::qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values = q.values_rdadapt))
}

#' Genotype-environement correlation test
#'
#' @param gen dosage matrix
#' @param var dataframe with predictor variables
#'
#' @return dataframe with r and p-values from correlation test
#' @export
rda_cor <- function(gen, var){
  cor_df <- purrr::map_dfr(colnames(gen), rda_cor_env_helper, gen, var)
  rownames(cor_df) <- NULL
  colnames(cor_df) <- c("r", "p", "snp", "var")
  return(cor_df)
}

#' Helper function for rda_cor_test
#'
#' @export
#' @noRd
rda_cor_env_helper <- function(snp_name, snp_df, env){
  cor_df <- data.frame(t(apply(env, 2, rda_cor_helper, snp_df[,snp_name])))
  cor_df$snp <- snp_name
  cor_df$env <- colnames(env)
  return(cor_df)
}

#' Helper function for rda_cor_test
#'
#' @export
#' @noRd
rda_cor_helper <- function(envvar, snp){
  if(sum(!is.na(envvar)) < 3 | sum(!is.na(snp)) < 3) return(c(r = NA, p = NA))
  mod <- stats::cor.test(envvar, snp, alternative = "two.sided", method = "pearson", na.action = "na.omit")
  pvalue <- mod$p.value
  r <- mod$estimate
  results <- c(r, pvalue)
  names(results) <- c("r", "p")
  return(results)
}

#' Plot RDA results
#' @param mod model object of class `rda`
#' @param rda_snps vector of outlier SNPs
#' @param pvalues if creating a manhattan plot (i.e., `manhattan = TRUE`), a matrix of p-values
#' @param axes which RDA axes to include while plotting (defaults to `all`)
#' @param biplot_axes if creating an RDA biplot (i.e., `rdaplot = TRUE`), which pairs of axes to plot. Defaults to plotting all pairs of axes possible, otherwise can be set to a single pair of axes (e.g., c(1,2)) or a list of axes pairs (e.g., list(c(1,2), c(2,3))))
#' @param manhattan whether to produce manhattan plot (defaults to `TRUE`)
#' @param rdaplot whether to produce an RDA biplot (defaluts to `TRUE`). If only one axes is provided, instead of a biplot a histogram will be created
#'
#' @export
#'
rda_plot <- function(mod, rda_snps, pvalues = NULL, axes = "all", biplot_axes = NULL, sig = 0.05, manhattan = TRUE, rdaplot = TRUE, binwidth = NULL){
  # Get axes
  if(axes == "all") axes <- 1:ncol(mod$CCA$v)

  # Make and get tidy dataframes for plotting
  tidy_list <- rda_ggtidy(mod, rda_snps, axes = axes)
  TAB_snps <- tidy_list[["TAB_snps"]]
  TAB_var <- tidy_list[["TAB_var"]]

  # Make RDA plots
  if(rdaplot){
    if(length(axes) == 1) {
      print(rda_hist(TAB_snps, binwidth = binwidth))
    } else if(!is.null(biplot_axes)){
      if(is.vector(biplot_axes)) print(rda_biplot(TAB_snps, TAB_var, biplot_axes = biplot_axes))
      if(is.list(biplot_axes)) lapply(biplot_axes, function(x) {print(rda_biplot(TAB_snps, TAB_var, biplot_axes = x))})
    } else {
      cb <- combn(length(axes), 2)
      if(!is.null(dim(cb))) {
        apply(cb, 2, function(x) {print(rda_biplot(TAB_snps, TAB_var, biplot_axes = x))})
      } else print(rda_biplot(TAB_snps, TAB_var, biplot_axes = cb))
    }
  }

  # Make manhattan plot
  if(manhattan & !is.null(pvalues)) print(rda_manhattan(TAB_snps, rda_snps, pvalues, sig = sig))
}


#' Make dataframe for ggplot from RDA results
#'
#' @export
#' @noRd
rda_ggtidy <- function(mod, rda_snps, axes){
  snp_scores <- vegan::scores(mod, choices = axes, display = "species", scaling = "none") # vegan references "species", here these are the snps
  TAB_snps <- data.frame(names = row.names(snp_scores), snp_scores)

  TAB_snps$type <- "Neutral"
  TAB_snps$type[TAB_snps$names %in% rda_snps] <- "Outliers"
  TAB_snps$type <- factor(TAB_snps$type, levels = c("Neutral", "Outliers"))
  TAB_var <- as.data.frame(vegan::scores(mod, choices = axes, display="bp")) # pull the biplot scores

  tidy_list <- list(TAB_snps = TAB_snps, TAB_var = TAB_var)
  return(tidy_list)
}


#' Helper function to plot RDA biplot
#'
#' @export
#' @noRd
rda_biplot <- function(TAB_snps, TAB_var, biplot_axes = c(1,2)){

  # Select axes for plotting
  xax <- paste0("RDA",biplot_axes[1])
  yax <- paste0("RDA",biplot_axes[2])
  TAB_snps_sub <- TAB_snps[,c(xax, yax, "type")]
  colnames(TAB_snps_sub) <- c("x", "y", "type")
  TAB_var_sub <- TAB_var[,c(xax, yax)]
  colnames(TAB_var_sub) <- c("x", "y")

  # scale the variable loadings for the arrows
  TAB_var_sub$x <- TAB_var_sub$x * max(TAB_snps_sub$x)/stats::quantile(TAB_var_sub$x)[4]
  TAB_var_sub$y <- TAB_var_sub$y * max(TAB_snps_sub$y)/stats::quantile(TAB_var_sub$y)[4]

  ## Biplot of RDA snps and variables scores
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    ggplot2::geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    ggplot2::geom_point(data = TAB_snps_sub, ggplot2::aes(x=x, y=y, colour = type), size = 1.4) +
    ggplot2::scale_color_manual(values = c(rgb(0.7,0.7,0.7,0.1), "#F9A242FF")) +
    ggplot2::geom_segment(data = TAB_var_sub, ggplot2::aes(xend=x, yend=y, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
    ggrepel::geom_text_repel(data = TAB_var_sub, ggplot2::aes(x=x, y=y, label = row.names(TAB_var_sub)), size = 4) +
    ggplot2::xlab(xax) + ggplot2::ylab(yax) +
    ggplot2::guides(color = ggplot2::guide_legend(title="snp type")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=ggplot2::rel(.8)),
                   strip.text = ggplot2::element_text(size=11))

}

#' Helper function to plot RDA manhattan plot
#'
#' @export
#' @noRd
rda_manhattan <- function(TAB_snps, rda_snps, pvalues, sig = 0.05){

  ## Manhattan plot
  TAB_manhattan <- data.frame(pos = 1:nrow(TAB_snps),
                              pvalues = pvalues,
                              type = factor(TAB_snps$type, levels = c("Neutral", "Outliers")))
  TAB_manhattan <- TAB_manhattan[order(TAB_manhattan$pos),]

  ggplot2::ggplot(data = TAB_manhattan) +
    ggplot2::geom_point(ggplot2::aes(x=pos, y=-log10(pvalues), col = type), size=1.4) +
    ggplot2:: scale_color_manual(values = c(rgb(0.7,0.7,0.7,0.5), "#F9A242FF", "#6B4596FF")) +
    ggplot2::xlab("position") + ggplot2::ylab("-log10(p)") +
    ggplot2::geom_hline(yintercept=-log10(sig), linetype="dashed", color = "black", size=0.6) +
    ggplot2::guides(color=ggplot2::guide_legend(title="snp type")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position="right",
                   legend.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=ggplot2::rel(.8)),
                   strip.text = ggplot2::element_text(size=11))
}

#' Helper function to plot RDA histogram
#'
#' @export
#' @noRd
rda_hist <- function(TAB_snps, binwidth = NULL){
  ggplot2::ggplot() +
    ggplot2::geom_histogram(data = TAB_snps, ggplot2::aes(fill = type, x = get(colnames(TAB_snps)[2])), binwidth = binwidth) +
    ggplot2::scale_fill_manual(values = c(rgb(0.7, 0.7, 0.7, 0.5), "#F9A242FF")) +
    ggplot2::guides(fill = ggplot2::guide_legend(title="snp type")) +
    ggplot2::xlab(colnames(TAB_snps)[2]) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="right",
                   legend.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=ggplot2::rel(.8)),
                   strip.text = ggplot2::element_text(size=11))
}

#' Make publication quality table of correlation test results
#'
#' @param cor_df dataframe of correlation results output from \link[algatr]{rda_cor}
#'
#' @return An object of class `gt_tbl`
#' @export
#'
rda_table <- function(cor_df, sig = 0.05, sig_only = TRUE, top = FALSE, order = FALSE, var = NULL, nrow = NULL, digits = 2){

  if(!is.null(var)) cor_df <- cor_df[cor_df$var %in% var, ]
  if(sig_only) cor_df <- cor_df[cor_df$p < sig, ]
  if(order) cor_df <- cor_df[order(abs(cor_df$r), decreasing = TRUE),]
  if(top) cor_df <- cor_df %>%
      dplyr::group_by(snp) %>%
      dplyr::filter(abs(r) == max(abs(r)))
  if(!is.null(nrow)) cor_df <- cor_df[1:nrow, ]

  cor_df <- cor_df %>% dplyr::as_tibble()
  if(!is.null(digits)) cor_df <- cor_df %>% dplyr::mutate(dplyr::across(-c(var, snp), round, digits))

  d <- max(abs(min(cor_df$r)), abs(max(cor_df$r)))

  suppressWarnings(
    tbl <- cor_df  %>%
      gt::gt() %>%
      gtExtras::gt_hulk_col_numeric(r, trim = TRUE, domain = c(-d,d))

  )

  tbl
}



