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
#  # For environmental and geographic data processing:
#  if (!require("RStoolbox", quietly = TRUE))devtools::install_github("bleutner/RStoolbox")
#  if (!require("geodata", quietly = TRUE)) install.packages("geodata")
#  if (!require("corrplot", quietly = TRUE)) install.packages("corrplot")
#  if (!require("vegan", quietly = TRUE)) install.packages("vegan")
#  if (!require("gdistance", quietly = TRUE)) install.packages("gdistance")
#  if (!require("topoDistance", quietly = TRUE)) install.packages("topoDistance")
#  if (!require("rmapshaper", quietly = TRUE)) install.packages("rmapshaper")
#  if (!require("wingen", quietly = TRUE)) devtools::github_install("wingen")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(terra)
library(raster)
library(RStoolbox)
library(ggplot2)
library(geodata)
library(viridis)
library(wingen)
library(tidyr)
library(tibble)

## ----test data----------------------------------------------------------------
load_algatr_example()

## ----get worldclim, eval=FALSE------------------------------------------------
#  wclim <- get_worldclim(coords = liz_coords, res = 5) # can set save_output = TRUE if you want a directory with files made

## ----wclim--------------------------------------------------------------------
# Load wclim object
data("wclim")
# Look at one of the layers, which corresponds to tile_15_wc2.1_30s_bio_1
wclim[[1]]

## ----wclim plot, warning = FALSE, fig.align='center', fig.height=5, fig.width=5----
plot(wclim[[1]], col = turbo(100), axes = FALSE)
points(liz_coords, pch = 19)

## ----collin env, results='hide'-----------------------------------------------
cors_env <- check_env(wclim)
# There are many pairs of layers that have correlation coefficients > 0.7 (the default threshold)

## ----collin env CA, fig.width = 10, fig.height = 10---------------------------
# For the sake of argument, we can also run this function on the CA_env object:
check_env(CA_env)
# None of the pairwise comparisons exceed the threshold value

## ----collin vals, fig.width = 10, fig.height = 10-----------------------------
# Once again, we can see that there are many pairs of variables that exceed the threshold
check_result <- check_vals(wclim, liz_coords)
head(check_result$cor_df)

## ----collin vals CA-----------------------------------------------------------
# Check using the CA_env object; note that these values are quite different from those from check_env()
check_vals(CA_env, liz_coords)

## ----collin dists-------------------------------------------------------------
check_results <- check_dists(wclim, liz_coords)
head(check_results$mantel_df)

## ----collin dists CA----------------------------------------------------------
# We can see that PC2 from CA_env and geo dists are significantly correlated!
check_dists(CA_env, liz_coords)

## ----raster PCA---------------------------------------------------------------
env_pcs <- rasterPCA(wclim, spca = TRUE)

# Let's take a look at the results for the top three PCs
plots <- lapply(1:3, function(x) ggR(env_pcs$map, x, geom_raster = TRUE))

plots[[1]] 
plots[[2]] 
plots[[3]]

## ----comp raster--------------------------------------------------------------
# We can also create a single composite raster plot with 3 PCs (each is assigned R, G, or B)
ggRGB(env_pcs$map, 1, 2, 3, stretch = "lin", q = 0)

## ----coords to raster---------------------------------------------------------
raster <- coords_to_raster(liz_coords, res = 1)
plot(raster, axes = FALSE)
points(liz_coords, pch = 19)

## ----extract enviro vars------------------------------------------------------
env <- raster::extract(CA_env, liz_coords)

# You can there are extracted values (n=53) for each environmental layer
head(env)

## ----env dist-----------------------------------------------------------------
env_dist <- env_dist(env)
plot(env_dist$CA_rPCA1)

## ----env dist heatmap---------------------------------------------------------
as.data.frame(env_dist$CA_rPCA1) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-"sample", names_to = "sample_comp", values_to = "dist") %>%
  ggplot(aes(x = as.numeric(sample), y = as.numeric(sample_comp), fill = dist)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis() +
  xlab("Sample") +
  ylab("Sample")

## ----geo dist-----------------------------------------------------------------
geo_dist <- geo_dist(liz_coords, type = "Euclidean")
plot(geo_dist)

## ----geo dist heatmap---------------------------------------------------------
# Make a fun heat map with the pairwise distances
geo_dist <- as.data.frame(geo_dist)
colnames(geo_dist) <- rownames(geo_dist)
geo_dist %>%
  rownames_to_column("sample") %>%
  gather("sample_comp", "dist", -"sample") %>%
  ggplot(aes(x = as.numeric(sample), y = as.numeric(sample_comp), fill = dist)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis() +
  xlab("Sample") +
  ylab("Sample")

## ----download DEM-------------------------------------------------------------
# Download the DEM raster using the elevation_30s function
# This will save a .tif file in your current dir
dem <- elevation_30s(country = "USA", path = getwd())

# Crop to California limits
CA_dem <- crop(dem, CA_env)

## ----plot dem-----------------------------------------------------------------
plot(CA_dem, axes = FALSE)
points(liz_coords, pch = 19)

## ----topo dists---------------------------------------------------------------
topo_dist <- geo_dist(liz_coords, type = "topo", lyr = CA_dem)

## ----topo dist heatmap--------------------------------------------------------
# Make a fun heat map with the pairwise distances
topo_dist <- as.data.frame(topo_dist)
colnames(topo_dist) <- 1:53
rownames(topo_dist) <- 1:53

topo_dist %>%
  rownames_to_column("sample") %>%
  gather("sample_comp", "dist", -"sample") %>%
  ggplot(aes(x = as.numeric(sample), y = as.numeric(sample_comp), fill = dist)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis() +
  xlab("Sample") +
  ylab("Sample")

## ----package vignette---------------------------------------------------------
vignette("topoDistance-vignette")
vignette("geodist")

