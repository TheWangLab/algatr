## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  # Install packages for installing other packages
#  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
#  
#  # For wingen:
#  if (!require("wingen", quietly = TRUE)) devtools::github_install("AnushaPB/wingen")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(wingen)
library(raster)
library(terra)

## ----test data----------------------------------------------------------------
load_algatr_example()
envlayer <- rast(CA_env$CA_rPCA1)

## ----coords to raster---------------------------------------------------------
liz_lyr <- coords_to_raster(liz_coords, res = 0.5, buffer = 5)

## ----preview gd, fig.width = 7, fig.height = 7, fig.align='center', warning=FALSE, message=FALSE----
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
preview_gd(liz_lyr,
  liz_coords,
  wdim = 3,
  fact = 0
)

## ----window gd, warning=FALSE-------------------------------------------------
wgd <- window_gd(liz_vcf,
  liz_coords,
  liz_lyr,
  stat = "pi",
  wdim = 3,
  fact = 0
)

## ----plot window, warning=FALSE, fig.width=7, fig.height=7, fig.align='center'----
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot map of pi
plot_gd(wgd, main = "Moving window pi", bkg = envlayer)
# Plot sample count map
plot_count(wgd, main = "Sample counts")

## ----krige window, warning=FALSE----------------------------------------------
kgd <- krig_gd(wgd, index = 1:2, liz_lyr, disagg_grd = 5)
summary(kgd)

## ----plot krige, fig.height=7, fig.width=7, fig.align='center', warning=FALSE----
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot kriged map of pi
plot_gd(kgd, main = "Kriged pi")
# Plot kriged sample count map
plot_count(kgd, main = "Kriged sample counts")

## ----mask sample count, fig.align='center', fig.width=7, fig.height=7, warning=FALSE----
mgd_1 <- mask_gd(kgd, kgd[["sample_count"]], minval = 1)
mgd_2 <- mask_gd(kgd, kgd[["sample_count"]], minval = 2)

par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
plot_gd(mgd_1, main = "Kriged & masked pi", bkg = envlayer)
plot_gd(mgd_2, main = "Kriged & masked pi", bkg = envlayer)

## ----mask sd, fig.height=7, fig.width=7, fig.align='center', warning=FALSE----
# Resample envlayer based on masked layer
r <- terra::resample(envlayer, mgd_1)

# Perform masking
mgd <- mask_gd(mgd_1, r)

# Plot masked map
plot_gd(mgd, main = "Kriged & masked pi", bkg = r)

## ----wingen do everything, warning=FALSE, results='hide'----------------------
results <- wingen_do_everything(
  gen = liz_vcf,
  lyr = liz_lyr,
  coords = liz_coords,
  wdim = 3,
  fact = 0,
  sample_count = TRUE,
  preview = FALSE,
  min_n = 2,
  stat = "pi",
  rarify = FALSE,
  kriged = TRUE,
  grd = liz_lyr,
  index = 1:2,
  agg_grd = NULL, disagg_grd = 4,
  agg_r = NULL, disagg_r = NULL,
  masked = TRUE, mask = envlayer,
  bkg = envlayer, plot_count = TRUE
)

## ----package vignette---------------------------------------------------------
# vignette("wingen-vignette")

