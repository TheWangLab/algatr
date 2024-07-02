# GDM with  custom matrices

# FORMAT DATA ---------------------------------------------------------------------------------------------------

load_algatr_example()
gendist <- liz_gendist
coords <- liz_coords

# Custom matrices (this should be a named list of matrices)
# e.g. list(custom1 = custom1_dist, custom2 = custom2_dist, etc.)
distPreds <- list(geodist = geo_dist(coords))

# Scale genetic distance data from 0 to 1
gendist <- scale01(gendist)

# Vector of sites (for individual-based sampling, this is just assigning 1 site to each individual)
site <- 1:nrow(gendist)

# Bind vector of sites with gen distances
gdmGen <- cbind(site, gendist)

# Bind vector of sites with distPreds
distPreds <- purrr::map(distPreds, ~as.matrix(cbind(site, .x)))

# Create dataframe of predictor variables
gdmPred <- data.frame(
  site = site,
  x = coords_df$x,
  y = coords_df$y
)

# If you have aedditional environmental data you can add it like this
env <- extract(CA_env, coords)
gdmPred <- data.frame(gdmPred, env)

# Format data for GDM
gdmData <- gdm::formatsitepair(gdmGen, bioFormat = 3, XColumn = "x", YColumn = "y", siteColumn = "site", predData = gdmPred)

# Add custom matrices
gdmData <- gdm::formatsitepair(gdmData, 4, predData = gdmPred, siteColumn = "site", distPreds = distPreds)
names(gdmData)

# By default, GDM just uses the index of the matrix in the list as the variable name (e.g. matrix_1 instead of geodist), which I find confusing, so I wrote this extra stuff below to rename the variables to your original names

# If the list of distPreds has unamed slots replace them 
distPreds_names <- names(distPreds)
if (any(distPreds_names == "")){
  empty_names <- distPreds_names[distPreds_names == ""]
  distPreds_names[distPreds_names == ""] <- paste0("matrix", seq_along(empty_names))
}

# Convert the gdmData names from matrix_1, matrix_2, etc. to the actual variable names
gdm_names <- paste0("matrix_", 1:length(distPreds_names))
key = data.frame(gdm_names = gdm_names, distPreds_names = distPreds_names)
# function to rename columns
rename_columns <- function(data, gdm_name, var_name) {
  names(data) <- ifelse(grepl(gdm_name, names(data)), sub(gdm_name, var_name, names(data)), names(data))
  return(data)
}
# Use purrr::reduce to apply it iteratively over the gdmData
gdmData <- purrr::reduce2(key$gdm_names, key$distPreds_names, rename_columns, .init = gdmData)
names(gdmData)

# RUN GDM -------------------------------------------------------------------------------------------------------

# Remove any remaining incomplete cases
cc <- stats::complete.cases(gdmData)
if (!all(cc)) {
  gdmData <- gdmData[cc, ]
  warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))
}

# Run GDM with all predictors
gdm_model_final <- gdm::gdm(gdmData, geo = TRUE)