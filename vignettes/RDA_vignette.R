## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  # Install packages for installing other packages
#  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#  
#  # For RDA:
#  if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
#  if (!require("qvalue", quietly = TRUE)) BiocManager::install("qvalue")
#  if (!require("robust", quietly = TRUE)) install.packages("robust")
#  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
#  if (!require("vegan", quietly = TRUE)) install.packages("vegan")

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(algatr)
library(dplyr)
library(raster)
library(vegan)

## ----gen data, warning = FALSE------------------------------------------------
load_algatr_example()
# Convert from vcf to dosage matrix:
gen <- vcf_to_dosage(liz_vcf)

## ----impute-------------------------------------------------------------------
# Are there NAs in the data?
gen[1:5, 1:5]
gen <- simple_impute(gen)
# Check that NAs are gone
gen[1:5, 1:5]

## ----enviro data--------------------------------------------------------------
# Extract environmental vars
env <- raster::extract(CA_env, liz_coords)

# Standardize environmental variables and make into dataframe
env <- scale(env, center = TRUE, scale = TRUE)
env <- data.frame(env)

## ----full model, warning = FALSE, message = FALSE-----------------------------
mod_full <- rda_run(gen, env, model = "full")

## ----full call----------------------------------------------------------------
mod_full$call

## ----full summary-------------------------------------------------------------
head(summary(mod_full))

## ----full Radj----------------------------------------------------------------
RsquareAdj(mod_full)

## ----best model, warning = FALSE, message = FALSE-----------------------------
mod_best <- rda_run(gen, env,
  model = "best",
  Pin = 0.05,
  R2permutations = 1000,
  R2scope = T
)

## ----best model stats---------------------------------------------------------
mod_best$call
mod_best$anova
RsquareAdj(mod_best)

## ----partial rda, warning=FALSE, message=FALSE--------------------------------
mod_pRDA_geo <- rda_run(gen, env, liz_coords,
  model = "full",
  correctGEO = TRUE,
  correctPC = FALSE
)

## ----pRDA results-------------------------------------------------------------
anova(mod_pRDA_geo)
RsquareAdj(mod_pRDA_geo) # 0.0305
head(summary(mod_pRDA_geo))

## ----pRDA struct+geo----------------------------------------------------------
mod_pRDA_gs <- rda_run(gen, env, liz_coords,
  model = "full",
  correctGEO = TRUE,
  correctPC = TRUE,
  nPC = 3
)

## ----pRDA struct+geo summaryl-------------------------------------------------
head(summary(mod_pRDA_gs))

## ----varpart, warning=FALSE, message=FALSE------------------------------------
varpart <- rda_varpart(gen, env, liz_coords,
  Pin = 0.05, R2permutations = 1000,
  R2scope = T, nPC = 3
)
rda_varpart_table(varpart, call_col = TRUE)

## ----pRDA final---------------------------------------------------------------
mod_pRDA <- rda_run(gen, env, model = "best", correctPC = TRUE, nPC = 2)

## ----pRDA summary-------------------------------------------------------------
mod_pRDA$anova

## ----loadings, warning=FALSE, message=FALSE-----------------------------------
rda_plot(mod_pRDA, axes = "all", binwidth = 20)

## ----z-scores-----------------------------------------------------------------
rda_sig_z <- rda_getoutliers(mod_pRDA, naxes = "all", outlier_method = "z", z = 3, plot = FALSE)

# How many outlier SNPs were detected?
length(rda_sig_z$rda_snps)

## ----cand snps----------------------------------------------------------------
rda_sig_p <- rda_getoutliers(mod_best, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.01, plot = FALSE)

# How many outlier SNPs were detected?
length(rda_sig_p$rda_snps)

## ----q-values-----------------------------------------------------------------
# Extract SNP names; choices is number of axes
snp_names <- rownames(scores(mod_best, choices = 2, display = "species"))

# Identify outliers that have q-values < 0.1
q_sig <-
  rda_sig_p$rdadapt %>%
  mutate(snp_names = snp_names) %>%
  filter(q.values <= 0.1)

# How many outlier SNPs were detected?
nrow(q_sig)

## ----intersect----------------------------------------------------------------
Reduce(intersect, list(
  q_sig$snp_names,
  rda_sig_p$rda_snps,
  rda_sig_z$rda_snps
))

## ----biplot-------------------------------------------------------------------
rda_plot(mod_best, rda_sig_p$rda_snps, biplot_axes = c(1, 2), rdaplot = TRUE, manhattan = FALSE)

# In the case of our partial RDA, there was only one RDA axis, so a histogram is generated
rda_plot(mod_pRDA, rda_sig_z$rda_snps, rdaplot = TRUE, manhattan = FALSE, binwidth = 0.01)

## ----Manhattan----------------------------------------------------------------
rda_plot(mod_best, rda_sig_p$rda_snps, rda_sig_p$pvalues, rdaplot = FALSE, manhattan = TRUE)

## ----simple results-----------------------------------------------------------
# Extract genotypes for outlier SNPs
rda_snps <- rda_sig_p$rda_snps
rda_gen <- gen[, rda_snps]

# Run correlation test
cor_df <- rda_cor(rda_gen, env)

# Make a table from these results (displaying only the first 5 rows):
rda_table(cor_df, nrow = 5)

# Order by the strength of the correlation
rda_table(cor_df, order = TRUE, nrow = 5)

# Only retain the top variable for each SNP based on the strength of the correlation
rda_table(cor_df, top = TRUE, nrow = 5)

# Display results for only one environmental variable
rda_table(cor_df, var = "CA_rPCA2", nrow = 5)

## ----simple RDA do everything, warning=FALSE, message=FALSE, results='asis', fig.align='center'----
results <- rda_do_everything(gen, CA_env, liz_coords,
  correctGEO = FALSE,
  correctPC = FALSE,
  outlier_method = "p",
  sig = 0.05,
  p_adj = "fdr",
  cortest = TRUE,
  varpart = FALSE
)

## ----pRDA do everything, warning=FALSE, message=FALSE, results='asis', fig.align='center'----
results <- rda_do_everything(gen, env, liz_coords,
  model = "full",
  correctGEO = TRUE,
  correctPC = TRUE,
  nPC = 3,
  varpart = TRUE,
  outlier_method = "z",
  z = 3,
  p_adj = "fdr"
)

