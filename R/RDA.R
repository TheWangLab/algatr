

#' Run RDA
#'
#' @param gen genotype matrix
#' @param env dataframe with environmental data
#' @param coords dataframe with coordinates (only needed if K selection is performed with TESS)
#' @param model whether to use all variables ("full") or perform variable selection ("best")
#' @param sig alpha level
#' @param padj_method correction method supplied to \code{p.adjust} (can be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#'
#' @inheritParams vegan::ordiR2step
#' @inheritParams stats::p.adjust
#' 
#' @return
#' @export
#'
#' @examples
rda_run <- function(gen, env, coords, model = "full", sig = 0.05, sig_method = "p", padj_method = "fdr", Pin = 0.05, R2permutations = 1000, R2scope = TRUE, lm_test = FALSE){
  
  # Modify environmental data --------------------------------------------------------------------------------------------------
  # Standardize environmental variables
  env <- scale(env, center = TRUE, scale = TRUE) 
  
  if(any(is.na(env))){
    warning("NA values found in env data, removing for RDA")
    gen <- gen[complete.cases(env),]
    env <- env[complete.cases(env),]
  }
  
  if(any(is.na(gen))){
    warning("NA values found in gen data, removing for RDA")
    gen <- gen[complete.cases(gen),]
    env <- gen[complete.cases(gen),]
  }
  
  # make into tibble
  env <- as_tibble(env)
  
  # Running RDA ----------------------------------------------------------------------------------------------------------------
  # Full model
  f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse='+')))
  RDAfull <- rda(f, data = env)
  
  if(model == "full"){
    final_mod <- RDAfull
  }
  
  if(model == "best"){
    # Null model
    RDA0 <- rda(gen ~ 1)
    
    # WEIRD ERROR FIX LATER
    CA_rPCA1 <- env$CA_rPCA1
    CA_rPCA2 <- env$CA_rPCA2
    CA_rPCA3 <- env$CA_rPCA3
    
    # Stepwise procedure with ordiR2step function
    final_mod <- vegan::ordiR2step(RDA0, RDAfull, Pin = Pin, R2permutations = R2permutations, R2scope = R2scope, trace = FALSE)
    
    # WRITE WARNING FOR IF BEST MODEL IS RDA0
    
  }

  # Identify candidate loci ----------------------------------------------------------------------------------------------------
  
  # Running the function with all axes
  # NOTE MODIFY THIS TO TAKE DIFFERENT NUMBERS OF AXES (USER DEFINED/AUTOMATICALLY CHOSEN USING SCREE or VAR explained)
  naxes <- ncol(final_mod$CCA$v)
  
  if(sig_method == "sd"){
    rda_loci_snp <- rda_sd(final_mod, K = naxes, z = z)$snp
    rda_loci <- which(colnames(gen) %in% rda_loci_snp)
  }
  
  if(sig_method == "p"){
    rdadapt_env <- rdadapt(final_mod, naxes)
    
    # P-values threshold after FDR correction (different from Capblancq & Forester 2021)
    pvalues <- p.adjust(rdadapt_env$p.values, method = padj_method)
    
    # Capblancq include a step where they only take pvalues with highest loading for each contig to deal with LD (not applied here)
    # NOTE: I Think this filtering step should occur before (e.g. only one loci per LD block, but you know which comes from where)
    
    ## Identifying the loci that are below the p-value threshold
    ##ADD THIS: . Selected loci can then be tested for their association with proposed environmental drivers using simple correlations  with  allele  frequencies,  or  permutations  (e.g.,  Pavlova  et al., 2013).
    #Identify rda cand loci (P)
    rda_loci <- which(pvalues < sig) 
  }
  
  # Identify environmental associations
  rda_gen <- gen[,rda_loci]
  if(lm_test){lm_df <- rda_lm_test(rda_gen, env)} else {lm_df <- NULL}
  
  # Plot results ---------------------------------------------------------------------------------------------------------------
  # Make and get tidy dataframes for plotting 
  tidy_list <- rda_tidy(final_mod, rda_loci, naxes)
  TAB_loci <- tidy_list[["TAB_loci"]]
  TAB_var <- tidy_list[["TAB_var"]]
  
  # Make RDA plots
  # ADJUST TO WORK FOR ALL AXES COMBOS
  print(rda_biplot(TAB_loci, TAB_var, rda_axes = c(1,2)))
  print(rda_biplot(TAB_loci, TAB_var, rda_axes = c(2,3)))
  
  # Make manhattan plot
  if(sig_method == "p"){print(rda_manhattan(TAB_loci, pvalues, sig = sig))}
  
  # Compile results ------------------------------------------------------------------------------------------------------------
  results <- list(rda_loci = rda_loci,
                  lm_df = lm_df,
                  rda_mod = final_mod)
  
  return(results)
}

# Function to conduct a RDA based genome scan from Capblancq & Forester 2021
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd
# NOTE: GO THROUGH THIS CODE
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

rda_sd <- function(rda, K, z = 3){
  #load scores
  load.rda <- scores(rda, choices=c(1:K), display="species")
  
  #identify loadings
  cand <- map_dfr(.x = 1:K, .f = cand_id, load.rda = load.rda, z = 3)
  colnames(cand) <- colnames(cand) <- c("axis","snp","loading")
  
  return(cand)
}

#OUTLIER FUNCTION
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#IDENTIFY CAND
cand_id <- function(i, load.rda, z = 3){
  cand <- outliers(load.rda[,i], z = z)
  cand <- bind_cols(rep(i, times=length(cand)), names(cand), unname(cand))
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
  env_p <- map_dfc(rda_gen, rda_lm_env_helper, env)
  rda_df <- bind_cols(colnames(env), env_p)
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
  env_p <- map_df(env, rda_lm_helper, locus = locus)
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
  mod_df <- bind_cols(locus, envvar)
  colnames(mod_df) <- c("locus", "envvar")
  mod <- lm(locus ~ envvar, data = mod_df)
  mod_sum <- summary(mod)
  pvalue <- mod_sum$coefficients[,"Pr(>|t|)"]["envvar"]
  return(pvalue)
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
rda_tidy <- function(rda_mod, rda_loci, naxes){
  locus_scores <- scores(rda_mod, choices = 1:naxes, display="species", scaling="none") # vegan references "species", here these are the loci
  TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
  TAB_loci$type <- "Neutral"
  TAB_loci$type[rda_loci] <- "Outliers"
  TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "Outliers"))
  TAB_loci <- TAB_loci[order(TAB_loci$type),]
  TAB_var <- as.data.frame(scores(rda_mod, choices=1:naxes, display="bp")) # pull the biplot scores
  
  tidy_list <- list(TAB_loci = TAB_loci, TAB_var = TAB_var)
  return(tidy_list)
}


#' Plot RDA biplot
#'
#' @param TAB_loci 
#' @param TAB_var 
#' @param rda_axes 
#'
#' @return
#' @export
#'
#' @examples
rda_biplot <- function(TAB_loci, TAB_var, rda_axes = c(1,2)){
  
  # Select axes for plotting
  xax <- paste0("RDA",rda_axes[1])
  yax <- paste0("RDA",rda_axes[2])
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
#' @param sig 
#'
#' @return
#' @export
#'
#' @examples
rda_manhattan <- function(TAB_loci, pvalues, sig){
  # get candidate loci
  rda_loci <- which(pvalues < sig) 
  
  ## Manhattan plot
  Outliers <- rep("Neutral", nrow(TAB_loci))
  Outliers[rda_loci] <- "Outliers"
  Outliers <- factor(Outliers, levels = c("Neutral", "Outliers"))
  TAB_manhatan <- data.frame(pos = 1:nrow(TAB_loci), 
                             pvalues = pvalues, 
                             Outliers = Outliers)
  TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
  ggplot(data = TAB_manhatan) +
    geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
    scale_color_manual(values = c(rgb(0.7,0.7,0.7,0.5), "#F9A242FF", "#6B4596FF")) +
    xlab("Loci") + ylab("-log10(p.values)") +
    geom_hline(yintercept=-log10(sig), linetype="dashed", color = "black", size=0.6) +
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
  tbl <-lm_df %>%
    mutate(across(-env, round, 3)) %>%
    gt() %>%
    data_color(columns = -env, colors = scales::col_numeric(
      palette = c("#F9A242FF"),
      domain = c(0, 0.05)
    )
    )
  tbl
}

#TODO: update docs and combine plot functions
