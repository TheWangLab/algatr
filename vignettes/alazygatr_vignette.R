## ----color, echo = FALSE, results = FALSE-------------------------------------
options(crayon.enabled = TRUE)

old_hooks <- fansi::set_knit_hooks(knitr::knit_hooks,
  which = c("output", "message", "error")
)

## ---- eval = FALSE------------------------------------------------------------
#  # Install packages for installing other packages
#  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
#  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#  
#  # For genetic distance processing
#  if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
#  if (!require("Assoctests", quietly = TRUE)) install.packages("Assoctests")
#  if (!require("readr", quietly = TRUE)) install.packages("readr")
#  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
#  if (!require("ecodist", quietly = TRUE)) install.packages("ecodist")
#  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
#  
#  # For genetic data processing
#  if (!require("gdsfmt", quietly = TRUE)) BiocManager::install("gdsfmt")
#  if (!require("SeqArray", quietly = TRUE)) BiocManager::install("SeqArray")
#  if (!require("SNPRelate", quietly = TRUE)) BiocManager::install("SNPRelate")
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
#  
#  # For LFMM:
#  if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
#  if (!require("Assoctests", quietly = TRUE)) install.packages("Assoctests")
#  if (!require("lfmm", quietly = TRUE)) install.packages("lfmm")
#  if (!require("TESS3_encho_sen", quietly = TRUE)) devtools::install_github("bcm-uga/TESS3_encho_sen")
#  if (!require("LEA", quietly = TRUE)) BiocManager::install("LEA")
#  if (!require("tess3r", quietly = TRUE)) install.packages("tess3r")
#  
#  # For RDA:
#  if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
#  if (!require("qvalue", quietly = TRUE)) BiocManager::install("qvalue")
#  if (!require("robust", quietly = TRUE)) install.packages("robust")
#  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
#  if (!require("vegan", quietly = TRUE)) install.packages("vegan")
#  
#  # For GDM:
#  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
#  if (!require("gdm", quietly = TRUE)) install.packages("gdm")
#  
#  # For MMRR:
#  if (!require("GGally", quietly = TRUE)) install.packages("GGally")
#  
#  # For TESS:
#  if (!require("automap", quietly = TRUE)) install.packages("automap")
#  if (!require("graphics", quietly = TRUE)) install.packages("graphics")
#  if (!require("LEA", quietly = TRUE)) BiocManager::install("LEA") # required by tess3r
#  if (!require("TESS3_encho_sen", quietly = TRUE)) devtools::install_github("bcm-uga/TESS3_encho_sen")
#  if (!require("LEA", quietly = TRUE)) BiocManager::install("LEA")
#  if (!require("fields", quietly = TRUE)) install.packages("fields")
#  if (!require("rworldmap", quietly = TRUE)) install.packages("rworldmap")
#  if (!require("tess3r", quietly = TRUE)) install.packages("tess3r")
#  
#  # For wingen:
#  if (!require("wingen", quietly = TRUE)) devtools::github_install("AnushaPB/wingen")
#  

## ----setup, include = FALSE, warning = FALSE, message = FALSE-----------------
library(algatr)

## ----load data----------------------------------------------------------------
load_algatr_example()
gen <- liz_vcf[, 1:21]
coords <- liz_coords[1:20, ]
envlayers <- CA_env

## ----alazygatr, warning = FALSE-----------------------------------------------
lazy_results <- do_everything_for_me(gen, coords, envlayers, quiet = FALSE)

