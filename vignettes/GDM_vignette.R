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
#  # For GDM:
#  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
#  if (!require("gdm", quietly = TRUE)) install.packages("gdm")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(raster)
library(terra)

## ----data---------------------------------------------------------------------
load_algatr_example()

## ----extract enviro vars------------------------------------------------------
env <- raster::extract(CA_env, liz_coords)

## ----gdm full-----------------------------------------------------------------
gdm_full <- gdm_run(
  gendist = liz_gendist,
  coords = liz_coords,
  env = env,
  model = "full",
  scale_gendist = TRUE
)

## ----full summary, warning = FALSE--------------------------------------------
summary(gdm_full$model)

## ----exp and obs dists, fig.align='center', fig.height=4, fig.width=8---------
gdm_plot_diss(gdm_full$model)

## ----isplines, fig.width = 3, fig.height = 3----------------------------------
par(mfrow = c(2, 2))
gdm_plot_isplines(gdm_full$model)

## ----GDM table----------------------------------------------------------------
gdm_table(gdm_full)

## ----GDM plots, results='hide'------------------------------------------------
gdm_map(gdm_full$model, CA_env, liz_coords)

## ----gdm masked, warning=FALSE, fig.align='center', fig.width=5, fig.height=5, warning = FALSE----
# Extract the GDM map from the GDM model object
map <- gdm_map(gdm_full$model, CA_env, liz_coords, plot_vars = FALSE)
maprgb <- map$pcaRastRGB

# Now, use `extrap_mask()` to do buffer-based masking
map_mask <- extrap_mask(liz_coords, maprgb, method = "buffer", buffer_width = 1.25)

# Plot the resulting masked map
plot_extrap_mask(maprgb, map_mask, RGB = TRUE, mask_col = rgb(1, 1, 1, alpha = 0.6))

## ----do everything, fig.width=5, fig.height=5, fig.align='center', results='hide'----
gdm_full_everything <- gdm_do_everything(liz_gendist,
  liz_coords,
  envlayers = CA_env,
  model = "full",
  scale_gendist = TRUE
)

