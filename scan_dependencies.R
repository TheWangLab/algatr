library(tidyr)
library(dplyr)
d <- renv::dependencies()
unique(d$Source)

methods <- c("LFMM", "RDA", "GDM", "MMRR", "TESS", "wingen", "masks", "get_worldclim", "Gen_dist", "SNP_LD", "check_vars", "dist", "format_data")
methods_d <- purrr::map(methods, \(x) {
  return(d %>% filter(grepl(x, Source)) %>% select(Package) %>% distinct() %>% mutate(method = x))
}) %>% bind_rows()

general_packages <-
  d %>%
  select(Package) %>%
  group_by(Package) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n > 5) %>%
  select(Package)

methods_d <- methods_d %>% filter(!(Package %in% general_packages$Package))

leftovers <-
  d %>%
  filter(!(Package %in% c(general_packages$Package, methods_d$Package))) %>%
  distinct(Package)

# tigris - used only in DATASET
# fansi - used only in alazygatr
# usethis, crayon, devtools, renv, BiocManager can all be SUGGESTS
# RStoolbox only used in the enviro_data_vignette
