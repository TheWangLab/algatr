## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  # Install packages for installing other packages
#  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
#  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#  
#  # For LFMM:
#  if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
#  if (!require("Assoctests", quietly = TRUE)) install.packages("Assoctests")
#  if (!require("lfmm", quietly = TRUE)) install.packages("lfmm")
#  if (!require("TESS3_encho_sen", quietly = TRUE)) devtools::install_github("bcm-uga/TESS3_encho_sen")
#  if (!require("LEA", quietly = TRUE)) BiocManager::install("LEA")
#  if (!require("tess3r", quietly = TRUE)) install.packages("tess3r")
#  

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(here)
library(raster)

## ----data, warning = FALSE----------------------------------------------------
load_algatr_example()
# Convert vcf to dosage matrix
gen <- vcf_to_dosage(liz_vcf)
# Also, our code assumes that sample IDs from gendist and coords are the same order; be sure to check this before moving forward!
env <- raster::extract(CA_env, liz_coords)

## ----impute-------------------------------------------------------------------
# Are there NAs in the data?
gen[1:5, 1:5]
gen <- simple_impute(gen)
# Check that NAs are gone
gen[1:5, 1:5]

## ----K selection, warning = FALSE, fig.align='center', fig.height=3, fig.width=3----
# Keep relevant params but retain default values for them
select_K(gen, K_selection = "tracy_widom", criticalpoint = 2.0234) # 6

select_K(gen, K_selection = "quick_elbow", low = 0.08, max.pc = 0.90) # 3

select_K(gen, K_selection = "tess", coords = liz_coords, Kvals = 1:10) # 3

select_K(gen, K_selection = "find_clusters", perc.pca = 90, max.n.clust = 10) # 4

## ----lfmm run, warning = FALSE, message = FALSE-------------------------------
ridge_results <- lfmm_run(gen, env, K = 6, lfmm_method = "ridge")
lasso_results <- lfmm_run(gen, env, K = 6, lfmm_method = "lasso")

## ----cand SNPs table, fig.width=5, fig.height=5, fig.align='center'-----------
# Build tables for each of our LFMM runs, displaying only significant SNPs and ordering according to effect size (B)
lfmm_table(lasso_results$df, order = TRUE)
lfmm_table(ridge_results$df, order = TRUE)

## ----modify tables, fig.width=5, fig.height=5, fig.align='center'-------------
lfmm_table(lasso_results$df, order = TRUE, var = "CA_rPCA2")

# Be aware that if significant SNPs < nrow, function will return NULL object
# lfmm_table(lasso_results$df, sig_only = FALSE, order = FALSE, nrow=10)

# Similarly, the same will occur if you try to specify a variable that is not significantly associated with a SNP
# lfmm_table(lasso_results$df, sig_only = TRUE, var = "CA_rPCA1")

## ----qqplot, fig.width=5, fig.height=5, fig.align='center'--------------------
lfmm_qqplot(lasso_results$df)

## ----manhattan lasso, fig.width=5, fig.height=5, fig.align='center'-----------
# As displayed in our table from above, only six SNPs are visible on the lasso method plots as outliers:
lfmm_manhattanplot(lasso_results$df, sig = 0.05)

## ----manhattan ridge, fig.width=5, fig.height=5, fig.align='center'-----------
lfmm_manhattanplot(ridge_results$df, sig = 0.05)

## ----do everything, warning = FALSE,  results = FALSE, message = FALSE, warning = FALSE, fig.height=4, fig.width = 6, fig.align='center'----
results <- lfmm_do_everything(gen = gen, env = CA_env, coords = liz_coords, lfmm_method = "lasso", K_selection = "tracy_widom")

## ----package vignette---------------------------------------------------------
vignette("lfmm")

