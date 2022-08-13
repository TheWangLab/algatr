# Code to create example data -------------------------------------------------------------

# load and save vcf subset
liz_vcf <- vcfR::read.vcfR("inst/extdata/populations_r20.haplotypes.filtered_m70_randomSNP.vcf")
# subset out 1000 SNPs
liz_vcf <- liz_vcf[1:1000,]

# Get CA shapefile
# download states from tigris
states <- tigris::states(cb = TRUE)
# reproject into wgs84 to match coordinates
states <- sf::st_transform(states, 4326)
# convert to SPDF
states <- sf::as_Spatial(states)
# subset out CA
CA <- states[which(states$NAME == "California"), "STUSPS"]

# load and save coords
liz_coords <- read.table("inst/extdata/Scelop.coord")
# rename cols
colnames(liz_coords) <- c("x", "y")
# get sample IDs from vcf data
liz_coords$ID <- colnames(liz_vcf@gt)[-1]
# create spatial coordinates
sp::coordinates(liz_coords) <- ~x+y
# add CRS
raster::crs(liz_coords) <- raster::crs(CA)
# only include coordinates within CA
liz_coords <- liz_coords[CA,]
# get IDS of coordinates within CA (to use to subset VCF)
IDS <- liz_coords$ID
# create dataframe
liz_coords <- data.frame(liz_coords)
# only keep x and y
liz_coords <- liz_coords[,c("x","y")]
usethis::use_data(liz_coords, overwrite = TRUE)

# subset vcf to match coords
index <- colnames(liz_vcf@gt) %in% IDS
# first col is format col
index[1] <- TRUE
# subset vcf to match coords
liz_vcf <- liz_vcf[ , index]
usethis::use_data(liz_vcf, overwrite = TRUE)

# check IDs match (remember first col is format)
stopifnot(colnames(liz_vcf@gt)[-1] == IDS)


# load env data
CA_env <- raster::stack(list.files("inst/extdata/PC_layers/", full.names = TRUE))
usethis::use_data(CA_env, overwrite = TRUE)

# process plink genetic distances
plink_file = here("data", "liz_test.dist")
plink_id_file = here("data", "liz_test.dist.id")
gendist <- as.data.frame(readr::read_tsv(plink_file, col_names = FALSE))
plink_names <- readr::read_tsv(plink_id_file, col_names = FALSE) %>%
  dplyr::select(-`X1`) %>%
  as.matrix()
# Assign row and col names according to sampleID
rownames(gendist) <- plink_names
colnames(gendist) <- plink_names
usethis::use_data(gendist, overwrite = TRUE)
