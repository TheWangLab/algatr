
# Code to create example data -------------------------------------------------------------

# Raw Sceloporus data files from:
# Bouzid N, Archie JW, Anderson RA, Grummer JA, Leaché AD (2022)
# Evidence for ephemeral ring species formation during the diversification history of western fence lizards (Sceloporus occidentalis)
# Molecular Ecology, 31: 620-631. doi: https://doi.org/10.1111/mec.15836

# liz_vcf -----------------------------------------------------------------

# Load and save vcf subset
liz_vcf <- vcfR::read.vcfR("inst/extdata/populations_r20.haplotypes.filtered_m70_randomSNP.vcf")
# Subsample 1000 SNPs
liz_vcf <- liz_vcf[1:1000,]


# CA shapefile ------------------------------------------------------------

# Get CA shapefile
# Download states from tigris
states <- tigris::states(cb = TRUE)
# Reproject into WGS84 to match coordinates
states <- sf::st_transform(states, 4326)
# Convert to SPDF
states <- sf::as_Spatial(states)
# Subset out CA
CA <- states[which(states$NAME == "California"), "STUSPS"]


# liz_coords --------------------------------------------------------------

# Load and save coords
liz_coords <- read.table("inst/extdata/Scelop.coord")
# Rename cols
colnames(liz_coords) <- c("x", "y")
# Get sample IDs from vcf data
liz_coords$ID <- colnames(liz_vcf@gt)[-1]
# Create spatial coordinates
sp::coordinates(liz_coords) <- ~x+y
# Add CRS
raster::crs(liz_coords) <- raster::crs(CA)
# Only include coordinates within CA
liz_coords <- liz_coords[CA,]
# Get IDS of coordinates within CA (to use to subset VCF)
IDS <- liz_coords$ID
# Create dataframe
liz_coords <- data.frame(liz_coords)
# Only keep x and y
liz_coords <- liz_coords[,c("x","y")]
usethis::use_data(liz_coords, overwrite = TRUE)

# Subset vcf to match coords
index <- colnames(liz_vcf@gt) %in% IDS
# First col is format col
index[1] <- TRUE
# Subset vcf to match coords
liz_vcf <- liz_vcf[ , index]
usethis::use_data(liz_vcf, overwrite = TRUE)

# Check IDs match (remember first col is format)
stopifnot(colnames(liz_vcf@gt)[-1] == IDS)


# CA_env ------------------------------------------------------------------

# Load env data
CA_env <- raster::stack(list.files("inst/extdata/PC_layers/", full.names = TRUE))
raster::writeRaster(CA_env, "inst/extdata/CA_env.tif", overwrite = TRUE)
CA_env <- raster::readAll(CA_env)
usethis::use_data(CA_env, overwrite = TRUE)


dos <- dos[complete.cases(dos),]
prcomp(~., data.frame(dos))


# liz_gendist -------------------------------------------------------------

# Process plink genetic distances
liz_gendist <- as.data.frame(readr::read_tsv("inst/extdata/liz_test.dist", col_names = FALSE))
plink_names <- readr::read_tsv("inst/extdata/liz_test.dist.id", col_names = FALSE) %>%
  dplyr::select(-`X1`) %>%
  as.matrix()

# Assign row and col names according to sampleID
rownames(liz_gendist) <- plink_names
colnames(liz_gendist) <- plink_names
usethis::use_data(liz_gendist, overwrite = TRUE)
