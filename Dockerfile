# Dockerfile modified from: http://haines-lab.com/post/2022-01-23-automating-computational-reproducibility-with-r-using-renv-docker-and-github-actions/
# Start with geospatial v 4.3.1
FROM rocker/geospatial:4.3.1

# Relabel docker
LABEL org.opencontainers.image.description=""

# Set default user
ENV DEFAULT_USER=algatr

# Use renv version 0.17.3
ENV RENV_VERSION 0.17.3

# Install renv
RUN Rscript -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN Rscript -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# Create a directory named after our project directory
WORKDIR /algatr

# Copy the lockfile over to the Docker image
COPY renv.lock renv.lock

# Install all R packages specified in renv.lock
RUN Rscript -e 'renv::restore()'

# Install the algatr package from GitHub
RUN Rscript -e 'remotes::install_github("TheWangLab/algatr@dependencies", build_vignettes = TRUE)'
