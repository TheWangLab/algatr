
# safely knit vignettes (captures errors)
safe_render <- purrr::safely(rmarkdown::render)

# get vector of vignette paths
vignettes <- list.files(here::here("vignettes"), full.names = TRUE, pattern = "*.Rmd")

# knit all vignettes
results <- purrr::map(vignettes, ~safe_render(.x, "html_document"))
names(results) <- basename(vignettes)

# pull out errors
errors <- map(results, "error") %>% discard(is.null)
print(errors)


