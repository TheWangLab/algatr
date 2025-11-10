# Package index

## Genetic data processing

- [`data_processing_packages()`](https://thewanglab.github.io/algatr/reference/data_processing_packages.md)
  : Install data processing packages

- [`gen_to_geno()`](https://thewanglab.github.io/algatr/reference/gen_to_geno.md)
  : Convert dosage matrix or vcf to geno type object (N.B.: this only
  works for diploids!)

- [`geno_to_dosage()`](https://thewanglab.github.io/algatr/reference/geno_to_dosage.md)
  : Convert lfmm/geno matrix to dosage matrix (N.B.: this only works for
  diploids!)

- [`ld_prune()`](https://thewanglab.github.io/algatr/reference/ld_prune.md)
  :

  ld_prune prunes SNPs based on linkage disequilibrium using `SNPRelate`
  and `SeqArray` packages

- [`simple_impute()`](https://thewanglab.github.io/algatr/reference/simple_impute.md)
  : Impute NA values NOTE: use extreme caution when using this form of
  simplistic imputation. We mainly provide this code for creating test
  datasets and highly discourage its use in analyses.

- [`str_impute()`](https://thewanglab.github.io/algatr/reference/str_impute.md)
  :

  Imputation of missing values using population structure inferred with
  [`LEA::snmf`](https://rdrr.io/pkg/LEA/man/main_sNMF.html)

- [`vcf_to_dosage()`](https://thewanglab.github.io/algatr/reference/vcf_to_dosage.md)
  : Convert a vcf to a dosage matrix

## Environmental data processing

- [`envirodata_packages()`](https://thewanglab.github.io/algatr/reference/envirodata_packages.md)
  : Install environmental and geographic data processing packages
- [`check_dists()`](https://thewanglab.github.io/algatr/reference/check_dists.md)
  : Check geographic and environmental distances for collinearity
- [`check_env()`](https://thewanglab.github.io/algatr/reference/check_env.md)
  : Check environmental layers for collinearity
- [`check_vals()`](https://thewanglab.github.io/algatr/reference/check_vals.md)
  : Check extracted values for collinearity
- [`env_dist()`](https://thewanglab.github.io/algatr/reference/env_dist.md)
  : Calculate distance between environmental vars
- [`geo_dist()`](https://thewanglab.github.io/algatr/reference/geo_dist.md)
  : Calculate geographic distance between coordinates
- [`get_worldclim()`](https://thewanglab.github.io/algatr/reference/get_worldclim.md)
  : Download and merge WorldClim data for study area
- [`rm_islands()`](https://thewanglab.github.io/algatr/reference/rm_islands.md)
  : Remove islands from mapping

## Masking

- [`extrap_mask()`](https://thewanglab.github.io/algatr/reference/extrap_mask.md)
  [`range_mask()`](https://thewanglab.github.io/algatr/reference/extrap_mask.md)
  [`sd_mask()`](https://thewanglab.github.io/algatr/reference/extrap_mask.md)
  [`buffer_mask()`](https://thewanglab.github.io/algatr/reference/extrap_mask.md)
  [`chull_mask()`](https://thewanglab.github.io/algatr/reference/extrap_mask.md)
  : Create raster mask based on coordinates
- [`masking_packages()`](https://thewanglab.github.io/algatr/reference/masking_packages.md)
  : Install masking packages

## Genetic distances

- [`gen_dist()`](https://thewanglab.github.io/algatr/reference/gen_dist.md)
  : Calculate genetic distances
- [`gen_dist_corr()`](https://thewanglab.github.io/algatr/reference/gen_dist_corr.md)
  : Plot the relationship between two distance metrics
- [`gen_dist_hm()`](https://thewanglab.github.io/algatr/reference/gen_dist_hm.md)
  : Make heatmap of genetic distances
- [`gen_dist_packages()`](https://thewanglab.github.io/algatr/reference/gen_dist_packages.md)
  : Install genetic distance packages

## TESS

- [`tess_barplot()`](https://thewanglab.github.io/algatr/reference/tess_barplot.md)
  : Create TESS barplot
- [`tess_col_default()`](https://thewanglab.github.io/algatr/reference/tess_col_default.md)
  : Create default TESS color palette
- [`tess_do_everything()`](https://thewanglab.github.io/algatr/reference/tess_do_everything.md)
  : TESS function to do everything
- [`tess_ggbarplot()`](https://thewanglab.github.io/algatr/reference/tess_ggbarplot.md)
  : Create TESS barplot using ggplot2
- [`tess_ggplot()`](https://thewanglab.github.io/algatr/reference/tess_ggplot.md)
  : ggplot of TESS results
- [`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md)
  : Krige admixture values
- [`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md)
  : Test multiple K values
- [`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md)
  : Create a custom legend for TESS maps
- [`tess_packages()`](https://thewanglab.github.io/algatr/reference/tess_packages.md)
  : Install TESS packages
- [`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
  : Plot all kriged Q values for each K
- [`bestK()`](https://thewanglab.github.io/algatr/reference/bestK.md) :
  Best K Selection based on cross entropy
- [`geom_tess()`](https://thewanglab.github.io/algatr/reference/geom_tess.md)
  : Create geom of TESS results that can be added to a ggplot object

## MMRR

- [`mmrr_df()`](https://thewanglab.github.io/algatr/reference/mmrr_df.md)
  : Make nice dataframe from MMRR results

- [`mmrr_do_everything()`](https://thewanglab.github.io/algatr/reference/mmrr_do_everything.md)
  : MMRR function to do everything

- [`mmrr_packages()`](https://thewanglab.github.io/algatr/reference/mmrr_packages.md)
  : Install MMRR packages

- [`mmrr_plot()`](https://thewanglab.github.io/algatr/reference/mmrr_plot.md)
  : Plot MMRR results

- [`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md)
  : Run MMRR and return model object

- [`mmrr_table()`](https://thewanglab.github.io/algatr/reference/mmrr_table.md)
  :

  Create `gt` table of MMRR results

- [`mmrr_var_sel()`](https://thewanglab.github.io/algatr/reference/mmrr_var_sel.md)
  : mmrr_var_sel performs MMRR with backward elimination variable
  selection

- [`MMRR()`](https://thewanglab.github.io/algatr/reference/MMRR.md) :
  MMRR performs Multiple Matrix Regression with Randomization analysis

- [`unfold()`](https://thewanglab.github.io/algatr/reference/unfold.md)
  : unfold converts the lower diagonal elements of a matrix into a
  vector

## GDM

- [`gdm_coeffs()`](https://thewanglab.github.io/algatr/reference/gdm_coeffs.md)
  : Get coefficients for each predictor

- [`gdm_df()`](https://thewanglab.github.io/algatr/reference/gdm_df.md)
  : Create dataframe of GDM results

- [`gdm_do_everything()`](https://thewanglab.github.io/algatr/reference/gdm_do_everything.md)
  : GDM function to do everything (fit model, get coefficients, make and
  save raster)

- [`gdm_format()`](https://thewanglab.github.io/algatr/reference/gdm_format.md)
  : Format Data for Generalized Dissimilarity Modeling (GDM)

- [`gdm_map()`](https://thewanglab.github.io/algatr/reference/gdm_map.md)
  : Make map from model

- [`gdm_packages()`](https://thewanglab.github.io/algatr/reference/gdm_packages.md)
  : Install GDM packages

- [`gdm_plot_diss()`](https://thewanglab.github.io/algatr/reference/gdm_plot_diss.md)
  : Plot compositional dissimilarity spline plots

- [`gdm_plot_isplines()`](https://thewanglab.github.io/algatr/reference/gdm_plot_isplines.md)
  : Plot I-splines for each variable

- [`gdm_plot_vars()`](https://thewanglab.github.io/algatr/reference/gdm_plot_vars.md)
  : Create a PCA plot for GDM

- [`gdm_run()`](https://thewanglab.github.io/algatr/reference/gdm_run.md)
  : Run GDM and return model object

- [`gdm_table()`](https://thewanglab.github.io/algatr/reference/gdm_table.md)
  :

  Create `gt` table of GDM results

- [`gdm_var_sel()`](https://thewanglab.github.io/algatr/reference/gdm_var_sel.md)
  : Get best set of variables from a GDM model

- [`gdm_varimp_table()`](https://thewanglab.github.io/algatr/reference/gdm_varimp_table.md)
  : Generate a Variable Importance Table for GDM Models

- [`scale01()`](https://thewanglab.github.io/algatr/reference/scale01.md)
  : Scale genetic distances from 0 to 1

- [`scaleRGB()`](https://thewanglab.github.io/algatr/reference/scaleRGB.md)
  : Scale three layers of environmental data to R, G, and B for mapping

## RDA

- [`rda_cor()`](https://thewanglab.github.io/algatr/reference/rda_cor.md)
  : Genotype-environment correlation test

- [`rda_do_everything()`](https://thewanglab.github.io/algatr/reference/rda_do_everything.md)
  : RDA function to do everything

- [`rda_getoutliers()`](https://thewanglab.github.io/algatr/reference/rda_getoutliers.md)
  : Get significant outliers from RDA model

- [`rda_packages()`](https://thewanglab.github.io/algatr/reference/rda_packages.md)
  : Install RDA packages

- [`rda_plot()`](https://thewanglab.github.io/algatr/reference/rda_plot.md)
  : Plot RDA results

- [`rda_run()`](https://thewanglab.github.io/algatr/reference/rda_run.md)
  : Run RDA

- [`rda_table()`](https://thewanglab.github.io/algatr/reference/rda_table.md)
  :

  Create `gt` table of RDA results

- [`rda_varpart()`](https://thewanglab.github.io/algatr/reference/rda_varpart.md)
  : Partial RDA variance partitioning

- [`rda_varpart_table()`](https://thewanglab.github.io/algatr/reference/rda_varpart_table.md)
  :

  Create `gt` table with RDA variance partitioning results

## LFMM

- [`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md)
  : Convert LFMM results into a tidy dataframe for downstream processing

- [`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md)
  : LFMM function to do everything

- [`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md)
  : LFMM Manhattan Plot

- [`lfmm_packages()`](https://thewanglab.github.io/algatr/reference/lfmm_packages.md)
  : Install LFMM packages

- [`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md)
  : LFMM QQplot

- [`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md)
  : Run LFMM

- [`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md)
  :

  Create `gt` table of LFMM results

- [`select_K()`](https://thewanglab.github.io/algatr/reference/select_K.md)
  [`select_K_tw()`](https://thewanglab.github.io/algatr/reference/select_K.md)
  [`select_K_elbow()`](https://thewanglab.github.io/algatr/reference/select_K.md)
  [`select_K_tess()`](https://thewanglab.github.io/algatr/reference/select_K.md)
  [`select_K_fc()`](https://thewanglab.github.io/algatr/reference/select_K.md)
  : K selection

- [`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md)
  : Quickly choose an elbow for a PC

- [`tw()`](https://thewanglab.github.io/algatr/reference/tw.md) :
  Tracy–Widom test

## Wingen

- [`wingen_do_everything()`](https://thewanglab.github.io/algatr/reference/wingen_do_everything.md)
  : wingen function to do everything (preview and generate moving window
  maps, krige, and mask)
- [`wingen_packages()`](https://thewanglab.github.io/algatr/reference/wingen_packages.md)
  : Install wingen packages

## Alazygatr

- [`alazygatr_packages()`](https://thewanglab.github.io/algatr/reference/alazygatr_packages.md)
  : Install alazygatr packages

- [`do_everything_for_me()`](https://thewanglab.github.io/algatr/reference/do_everything_for_me.md)
  :

  Lazy run of all landscape genomic analyses contained within `algatr`

## Data

- [`load_algatr_example()`](https://thewanglab.github.io/algatr/reference/load_algatr_example.md)
  : Load example data
- [`CA_env`](https://thewanglab.github.io/algatr/reference/CA_env.md) :
  Example environmental data, calculated by performing a raster PCA on
  18 bioclimatic variables for state of California
- [`liz_coords`](https://thewanglab.github.io/algatr/reference/liz_coords.md)
  : Example coordinates from Bouzid et al. 2022
- [`liz_gendist`](https://thewanglab.github.io/algatr/reference/liz_gendist.md)
  : Example genetic distance matrix, calculated with Plink using data
  from Bouzid et al. 2022
- [`liz_vcf`](https://thewanglab.github.io/algatr/reference/liz_vcf.md)
  : Example VCF from Bouzid et al. 2022

## Other

- [`coords_to_sf()`](https://thewanglab.github.io/algatr/reference/coords_to_sf.md)
  : Convert from matrix, data frame, or sf to sf (sf is a pass through)
- [`coords_to_sp()`](https://thewanglab.github.io/algatr/reference/coords_to_sp.md)
  : Convert from matrix, data frame, or sf to formatted sp
