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
#  # For genetic data processing
#  if (!require("gdsfmt", quietly = TRUE)) BiocManager::install("gdsfmt")
#  if (!require("SeqArray", quietly = TRUE)) BiocManager::install("SeqArray")
#  if (!require("SNPRelate", quietly = TRUE)) BiocManager::install("SNPRelate")
#  if (!require("LEA", quietly = TRUE)) install.packages("LEA")
#  if (!require("ecodist", quietly = TRUE)) install.packages("ecodist")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(algatr)
library(here)
library(purrr)
library(dplyr)
library(ggplot2)

## ----vcf----------------------------------------------------------------------
# Load the example data
load_algatr_example()

# Look at genotypes for five individuals at eight sites:
liz_vcf@gt[1:8, 2:6]

## ----dosage-------------------------------------------------------------------
dosage <- vcf_to_dosage(liz_vcf)

# Look at genotypes for five individuals at five sites:
dosage[1:5, 1:8]

## ----simple impute------------------------------------------------------------
# Are there NAs in the data?
dosage[1:5, 1:5]
simple_dos <- simple_impute(dosage)
# Check that NAs are gone
simple_dos[1:5, 1:5]

## ----str impute, message = FALSE, results = FALSE, warning = FALSE------------
str_dos <- str_impute(gen = dosage, K = 1:3, entropy = TRUE, repetitions = 3, quiet = FALSE, save_output = FALSE)

## ----compare impute, message = FALSE, warning = FALSE, fig.align = 'center'----
str_dist <- as.matrix(ecodist::distance(str_dos, method = "euclidean"))
simple_dist <- as.matrix(ecodist::distance(simple_dos, method = "euclidean"))
gen_dist_corr(dist_x = simple_dist, dist_y = str_dist, metric_name_x = "simple", metric_name_y = "SNMF")

## ----package vignette---------------------------------------------------------
vignette("SNPRelate")

