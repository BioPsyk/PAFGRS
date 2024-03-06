# Use an official R base image
FROM rocker/r-ubuntu:22.04

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libpng-dev \
    libtiff5-dev \
    libudunits2-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Install R packages that are dependencies for your package and devtools
RUN R -e "install.packages('ggplot2')"  
RUN R -e "install.packages(c('data.table', 'Matrix', 'MASS', 'mvtnorm'))"  
RUN R -e "install.packages(c('methods', 'knitr', 'rmarkdown','testthat'))" 
RUN R -e "install.packages('devtools')" 
RUN R -e "install.packages('igraph')"
RUN R -e "install.packages('kinship2')"

# Set the working directory to a neutral location like /data for runtime
WORKDIR /pafgrs

RUN echo "library(devtools)" >> /usr/lib/R/etc/Rprofile.site
RUN echo "setwd('/pafgrs')" >> /usr/lib/R/etc/Rprofile.site
RUN echo "load_all()" >> /usr/lib/R/etc/Rprofile.site
ENTRYPOINT ["R", "--no-save"]


