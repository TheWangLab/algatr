# Start from the geospatial base image
FROM ghcr.io/rocker-org/devcontainer/geospatial:4.5

# Set default CRAN repo
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site

# Relax compiler warnings for R packages
ENV CFLAGS="-Wno-error=format-security"
ENV CXXFLAGS="-Wno-error=format-security"

# Install system dependencies for R packages (digest, shiny, etc.)
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libicu-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    pandoc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Avoid sandbox issues with renv or devtools
ENV RENV_CONFIG_SANDBOX_ENABLED=false
ENV R_COMPILE_AND_INSTALL_PACKAGES=always

# Install remotes for GitHub installation
RUN Rscript -e 'install.packages("remotes")'

# Install algatr and its dependencies
RUN Rscript -e 'remotes::install_github("TheWangLab/algatr", build_vignettes = FALSE)'
RUN Rscript -e 'algatr::alazygatr_packages()'
RUN Rscript -e 'remotes::install_github("TheWangLab/algatr", build_vignettes = TRUE)'
