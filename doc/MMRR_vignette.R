## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  # Install packages for installing other packages
#  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
#  
#  # For MMRR:
#  if (!require("GGally", quietly = TRUE)) install.packages("GGally")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(here)
library(raster)

## ----genetic data, warning = FALSE--------------------------------------------
load_algatr_example()
Y <- as.matrix(liz_gendist)

## ----enviro data--------------------------------------------------------------
# Extract enviro vars
env <- raster::extract(CA_env, liz_coords)
# Calculate environmental distances
X <- env_dist(env)
# Add geographic distance to X
X[["geodist"]] <- geo_dist(liz_coords)

## ----mmrr full----------------------------------------------------------------
set.seed(10)
results_full <- mmrr_run(Y, X, nperm = 99, stdz = TRUE, model = "full")

## ----mmrr best----------------------------------------------------------------
# Run MMRR with all variables
set.seed(01)
results_best <- mmrr_run(Y, X, nperm = 99, stdz = TRUE, model = "best")

## ----mmrr plots, fig.width = 5, fig.height = 5, fig.align='center'------------
# Single variable plot
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "vars", stdz = TRUE)

## ----fitted plot, fig.width = 5, fig.height = 5, fig.align='center'-----------
# Fitted variable plot
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "fitted", stdz = TRUE)

## ----cov plot, fig.width = 5, fig.height = 5, fig.align='center'--------------
# Covariance plot
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "cov", stdz = TRUE)

## ----mmrr best plots, fig.width = 5, fig.height = 5, fig.align='center'-------
mmrr_plot(Y, X, mod = results_best$mod, plot_type = "all", stdz = TRUE)

## ----stats--------------------------------------------------------------------
mmrr_table(results_full, digits = 2, summary_stats = TRUE)

## ----all vars, fig.width = 5, fig.height = 5, fig.align = 'center', warning = FALSE----
set.seed(01)
mmrr_full_everything <- mmrr_do_everything(liz_gendist, liz_coords, env = CA_env, geo = TRUE, model = "full")

## ----toy resist surface, fig.width = 5, fig.height = 5, fig.align = 'center', warning = FALSE----
# here we aggregate the layer for computational speed
lyr <- aggregate(CA_env$CA_rPCA1, 50)
plot(lyr)
points(liz_coords)

# Recreate MMRR input with resistance distances
# Calculate environmental distances
X <- env_dist(env)
# Add geographic distance to X
X[["resistdist"]] <- geo_dist(liz_coords, type = "resistance", lyr = lyr)

## ----mmrr resist, fig.width = 5, fig.height = 5, fig.align = 'center', warning = FALSE----
# Run MMRR with resistance distances
results_resist <- mmrr_run(Y, X, nperm = 99, stdz = TRUE, model = "full")
mmrr_plot(Y, X, mod = results_resist$mod, plot_type = "all", stdz = TRUE)
mmrr_table(results_resist)

