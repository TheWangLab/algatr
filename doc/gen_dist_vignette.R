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
#  # For genetic data processing (genetic distance only)
#  if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
#  if (!require("Assoctests", quietly = TRUE)) install.packages("Assoctests")
#  if (!require("readr", quietly = TRUE)) install.packages("readr")
#  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
#  if (!require("ecodist", quietly = TRUE)) install.packages("ecodist")
#  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(cowplot)

## ----load example-------------------------------------------------------------
load_algatr_example()

## ----euc dists----------------------------------------------------------------
# Calculate Euclidean distance matrix
euc_dists <- gen_dist(liz_vcf, dist_type = "euclidean")

euc_dists[1:5, 1:5]

## ----plink dists--------------------------------------------------------------
plink_dists <- gen_dist(plink_file = system.file("extdata", "liz_test.dist", package = "algatr"), plink_id_file = system.file("extdata", "liz_test.dist.id", package = "algatr"), dist_type = "plink")

## ----pc dists-----------------------------------------------------------------
pc_dists <- gen_dist(liz_vcf, dist_type = "pc", npc_selection = "auto", criticalpoint = 2.0234)

## ----gendist corr, fig.align='center', fig.height=5, fig.width=5--------------
# Plot some pairwise comparisons, providing names for the metrics
p_euc_plink <- gen_dist_corr(euc_dists, plink_dists, "Euclidean", "Plink")
p_pc_plink <- gen_dist_corr(pc_dists, plink_dists, "PC_based", "Plink")

# Show all plots as panels on a single plot
plot_grid(p_euc_plink, p_pc_plink, nrow = 1)

## ----euc hm, fig.align='center', fig.height=5, fig.width=5--------------------
gen_dist_hm(euc_dists)

## ----dps hm, fig.align='center', fig.height=5, fig.width=5--------------------
dps_dists <- gen_dist(liz_vcf, dist_type = "dps")
gen_dist_hm(dps_dists)

