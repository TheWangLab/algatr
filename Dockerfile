# Start from the geospatial base image
FROM ghcr.io/rocker-org/devcontainer/geospatial:4.5

# Default CRAN
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site

# Toolchain & headers commonly needed by CRAN packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libicu-dev libharfbuzz-dev libfribidi-dev \
    pandoc && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Relax compiler warnings (avoid -Werror issues)
ENV CFLAGS="-O2"
ENV CXXFLAGS="-O2"

# Avoid sandbox hiccups; always compile if needed
ENV RENV_CONFIG_SANDBOX_ENABLED=false
ENV R_COMPILE_AND_INSTALL_PACKAGES=always

# Preinstall remotes (optional convenience)
RUN Rscript -e 'install.packages("remotes", quiet = TRUE)'

# (Optional) Prewarm cache with your package/deps; skip vignettes first
RUN Rscript -e 'remotes::install_github("TheWangLab/algatr", build_vignettes = FALSE)'
RUN Rscript -e 'algatr::alazygatr_packages()'
RUN Rscript -e 'remotes::install_github("TheWangLab/algatr", build_vignettes = TRUE)"