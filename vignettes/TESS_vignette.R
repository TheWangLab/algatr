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
#  # For TESS:
#  if (!require("automap", quietly = TRUE)) install.packages("automap")
#  if (!require("graphics", quietly = TRUE)) install.packages("graphics")
#  if (!require("TESS3_encho_sen", quietly = TRUE)) devtools::install_github("bcm-uga/TESS3_encho_sen")
#  if (!require("LEA", quietly = TRUE)) BiocManager::install("LEA")
#  if (!require("fields", quietly = TRUE)) install.packages("fields")
#  if (!require("rworldmap", quietly = TRUE)) install.packages("rworldmap")
#  if (!require("tess3r", quietly = TRUE)) install.packages("tess3r")
#  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(here)
library(wingen)
library(tess3r)
library(ggplot2)
library(terra)
library(raster)
library(fields)
library(rworldmap)
library(automap)
library(cowplot)

# for plotting (from TESS)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

## ----data---------------------------------------------------------------------
load_algatr_example()
# Our code assumes that the first column is longitude and second is latitude; check this:
head(liz_coords)
# Also, our code assumes that sample IDs from gendist and coords are the same order; be sure to check this before moving forward!

## ----dosage-------------------------------------------------------------------
# Convert vcf to genotype matrix
liz_dosage <- vcf_to_dosage(liz_vcf)

## ----raster, fig.align = 'center'---------------------------------------------
# First, create a grid for kriging
# We can use one environmental layer (PC1), aggregated (i.e., increased cell size) to increase computational speed
krig_raster <- raster::aggregate(CA_env[[1]], fact = 6)

# If you want to see the difference between the non-aggregated (original) and aggregated rasters:
terra::plot(CA_env[[1]], col = mako(100), axes = FALSE)
terra::plot(krig_raster, col = mako(100), axes = FALSE)

## ----auto K selection, message = FALSE, fig.align = 'center'------------------
# Best K is 3; this provides a more reasonable estimate for the "best" K compared to manual selection above
tess3_result <- tess_ktest(liz_dosage, liz_coords, Kvals = 1:10, ploidy = 2, K_selection = "auto")

## ----no K selection, message = FALSE------------------------------------------
tess3_obj_noK <- tess3(liz_dosage, coord = as.matrix(liz_coords), K = 3, method = "projected.ls", ploidy = 2)

## ----qmatrix------------------------------------------------------------------
# Get TESS object and best K from results
tess3_obj <- tess3_result$tess3_obj
bestK <- tess3_result[["K"]]

# Get Qmatrix with ancestry coefficients
qmat <- qmatrix(tess3_obj, K = bestK)

# qmat contains ancestry coefficient values for each individual (row) and each K value (column)
head(qmat)

## ----barplot, fig.align = 'center'--------------------------------------------
tess_barplot(qmat)

## ----ggbarplot, fig.align = 'center'------------------------------------------
tess_ggbarplot(qmat, legend = FALSE)

