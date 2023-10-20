library(tidyr)
library(dplyr)
d <- renv::dependencies()
unique(d$Source)

# / in front of dist because otherwise Gen_dist packages get counted for dist
methods <- c("LFMM", "RDA", "GDM", "MMRR", "TESS", "wingen", "masks", "get_worldclim", "Gen_dist", "SNP_LD", "check_vars", "/dist", "format_data")
methods_d <- purrr::map(methods, \(x) {
  return(d %>% filter(grepl(x, Source)) %>% dplyr::select(Package) %>% distinct() %>% mutate(method = x))
}) %>% bind_rows()

general_packages <-
  d %>%
  dplyr::select(Package) %>%
  group_by(Package) %>%
  count() %>%
  filter(n > 4) %>%
  dplyr::select(Package) %>%
  # remove non-general packages
  filter(!(Package %in% c("wingen"))) %>%
  arrange(Package)

methods_d <- methods_d %>% filter(!(Package %in% general_packages$Package))
view(methods_d)

leftovers <-
  d %>%
  filter(!(Package %in% c("algatr", "wingen", "knitr", "rmarkdown"))) %>%
  filter(!(Package %in% c(general_packages$Package, methods_d$Package))) %>%
  distinct(Package)

print(leftovers)

# tigris - used only in DATASET
# fansi - used only in alazygatr
# usethis, crayon, devtools, renv, BiocManager - all SUGGESTS
