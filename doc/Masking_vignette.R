## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(algatr)
library(raster)
library(terra)
library(viridis)

## ----data, fig.align='center', fig.height=5, fig.width=5----------------------
# Load test data, including CA_env which are the envlayers we'll be using
load_algatr_example()

# For the purposes of simplicity, let's just use one of the PCs for mapping:
envlayers <- CA_env$CA_rPCA1

par(mar = c(0, 0, 0, 0))
# Let's take a look at the map with no masking:
plot(envlayers, col = viridis(100), axes = FALSE, box = FALSE)

## ----range mask, fig.align='center', fig.height=5, fig.width=5----------------
par(mar = c(0, 0, 0, 0))

# Extrapolate env values for given coordinates
map_mask <- extrap_mask(liz_coords, envlayers, method = "range")

# Now, plot the map with masked areas
plot_extrap_mask(envlayers, map_mask, RGB_cols = FALSE)

# Let's make the masked areas white with no transparency
plot_extrap_mask(envlayers, map_mask, RGB_cols = FALSE, mask_col = rgb(1, 1, 1, alpha = 1))

## ----sd mask, fig.align='center', fig.height=5, fig.width=5-------------------
par(mar = c(0, 0, 0, 0))

# Let's start with nsd=2
map_mask <- extrap_mask(liz_coords, envlayers, method = "sd", nsd = 2)
plot_extrap_mask(envlayers, map_mask, RGB = FALSE)

# Now, increase nsd to 3 and see how the map masking changes:
map_mask <- extrap_mask(liz_coords, envlayers, method = "sd", nsd = 3)
plot_extrap_mask(envlayers, map_mask, RGB = FALSE)

## ----buffer mask, fig.align='center', fig.height=5, fig.width=5---------------
par(mar = c(0, 0, 0, 0))

map_mask <- extrap_mask(liz_coords, envlayers, method = "buffer", buffer_width = 0.25)
plot_extrap_mask(envlayers, map_mask, RGB = FALSE)

# Increase buffer size
map_mask <- extrap_mask(liz_coords, envlayers, method = "buffer", buffer_width = 0.5)
plot_extrap_mask(envlayers, map_mask, RGB = FALSE)

# Increase buffer size and change masking color and transparency
map_mask <- extrap_mask(liz_coords, envlayers, method = "buffer", buffer_width = 1)
plot_extrap_mask(envlayers, map_mask, RGB = FALSE, mask_col = rgb(1, 1, 1, alpha = 1))

## ----chull mask, fig.align='center', fig.height=5, fig.width=5----------------
par(mar = c(0, 0, 0, 0))

map_mask <- extrap_mask(liz_coords, envlayers, method = "chull")
plot_extrap_mask(envlayers, map_mask, RGB = FALSE)

# Increase the buffer size
map_mask <- extrap_mask(liz_coords, envlayers, method = "chull", buffer_width = 0.5)
plot_extrap_mask(envlayers, map_mask, RGB = FALSE)

# Increase the buffer size again
map_mask <- extrap_mask(liz_coords, envlayers, method = "chull", buffer_width = 1)
plot_extrap_mask(envlayers, map_mask, RGB = FALSE)

## ----cali polygon, fig.align='center', fig.height=5, fig.width=5--------------
par(mar = c(0, 0, 0, 0))

states <- getData("GADM", country = "United States", level = 1)
# states <- gadm("United States", level = 1, path = here())
cali <- states[states$NAME_1 == "California", ]
plot(cali)

## ----remove islands, fig.align='center', fig.height=5, fig.width=5------------
par(mar = c(0, 0, 0, 0))

cali_noislands <- rm_islands(envlayers, cali)

# Now, let's plot the env layer and see how it's removed the Channel Islands
plot(cali_noislands)

## ----remove islands multilayers, fig.align='center', fig.height=5, fig.width=5----
par(mar = c(0, 0, 0, 0))

cali_noislands <- rm_islands(CA_env, cali)

# Islands are gone from all three enviro PC layers
plot(cali_noislands)

