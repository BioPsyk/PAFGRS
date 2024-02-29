# Use an official R base image
FROM r-base:4.1.0

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libudunits2-dev \
    libgdal-dev \
    libproj-dev \
    libcairo2-dev \
    libglpk-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages that are dependencies for your package and devtools
RUN R -e "install.packages(c('devtools', 'methods', 'knitr', 'rmarkdown','testthat'))" 
RUN R -e "install.packages(c('data.table', 'Matrix', 'MASS', 'mvtnorm'))"  
RUN R -e "install.packages('ggplot2')"  
RUN R -e "install.packages('igraph')"

# Copy the PAFGRS package files into the Docker image
COPY . /usr/local/src/PAFGRS

# Set the working directory to where we copied PAFGRS
WORKDIR /usr/local/src/PAFGRS

# Install the PAFGRS package using devtools
# Ensure dependencies are installed
RUN R -e "devtools::install_deps('/usr/local/src/PAFGRS')"
RUN R -e "devtools::install('/usr/local/src/PAFGRS', dependencies=TRUE)"

# Set the working directory to a neutral location like /data for runtime
WORKDIR /data

# Optionally, keep the container running if you're not specifying a command
CMD ["tail", "-f", "/dev/null"]

