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

## ---- cache = TRUE, results = FALSE, warning = FALSE--------------------------
krig_admix <- tess_krig(qmat, liz_coords, krig_raster)

## ----barplot, fig.align = 'center'--------------------------------------------
tess_barplot(qmat)

## ----ggbarplot, fig.align = 'center'------------------------------------------
tess_ggbarplot(qmat, legend = FALSE)

## ----basic ggplot, fig.width = 5, fig.height = 5, fig.align = 'center'--------
par(mfrow = c(2, 2), pty = "s", mar = rep(0, 4))

tess_ggplot(krig_admix, plot_method = "maxQ")
tess_ggplot(krig_admix, plot_method = "allQ", minQ = 0.20)
tess_ggplot(krig_admix, plot_method = "maxQ_poly")
tess_ggplot(krig_admix, plot_method = "allQ_poly", minQ = 0.20)

## ----extended mapping, fig.width = 5, fig.height = 5, warning = FALSE, fig.align = 'center'----
tess_ggplot(krig_admix, plot_method = "maxQ", ggplot_fill = scale_fill_manual(values = c("#bd9dac", "#257b94", "#476e9e")), plot_axes = TRUE)

## ----plot all K, fig.width = 10, fig.height = 4, fig.align = 'center', warning = FALSE----
par(mfrow = c(1, nlyr(krig_admix)), mar = rep(2, 4), oma = rep(1, 4))
tess_plot_allK(krig_admix, col_breaks = 20, legend.width = 2)

## ----tess default barplot, fig.align = 'center', message = FALSE--------------
barplot(qmat, sort.by.Q = TRUE, border = NA, space = 0, xlab = "Individuals", ylab = "Ancestry coefficients")

## ---- fig.width = 5, fig.height = 5, fig.align = 'center', message = FALSE, warning = FALSE----
plot(qmat,
  liz_coords,
  method = "map.max",
  interpol = FieldsKrigModel(10),
  main = "Ancestry coefficients",
  xlab = "Longitude", ylab = "Latitude",
  col.palette = CreatePalette(),
  resolution = c(300, 300), cex = .4
)

## ----do everything, cache = TRUE, results = FALSE, warning = FALSE, fig.width = 5, fig.height = 5, fig.align = 'center'----
# One could also use krig_raster as the grid object
results <- tess_do_everything(liz_vcf, liz_coords, grid = CA_env[[1]], Kvals = 1:10, K_selection = "auto")

