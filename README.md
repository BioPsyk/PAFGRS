# PA-FGRS

## Comprehensive Guide
An in-depth tutorial elucidating the usage and functionalities of PA-FGRS is available in the Articles section at https://biopsyk.github.io/PAFGRS

For those who wish to delve deeper and regenerate the vignette, the requisite code is provided below. Please execute this in your R environment for the desired output.

``` r
# Make sure to clone this repository
git clone git@github.com:BioPsyk/PAFGRS.git

# Enter the cloned repository
cd PAFGRS

# Start R and run
library(devtools)
build_vignettes()
```

The vignette in the format .html, .R and .RMD can now be found in the `doc` folder. 

## Containers
For the ability to run PAFGRS in places where new installations are not possible we provide a container image using docker and singularity. First use Docker to build a docker image on a machine where you have root access. Then use singularit to convert the docker image to a singularity image, which is a file that easily can be moved and copied. 

```
# Build docker image
docker build -t ibp-pafgrs .
# Convert to Singularity image
docker save ibp-pafgrs:latest -o ibp-pafgrs.tar
singularity build ibp-pafgrs.sif docker-archive://ibp-pafgrs.tar

# set folder to mount, e.g., your home folder
folder_to_mount=${HOME}
# Run the image starting R directly, mounted data accessible in /data inside the image.
singularity exec --bind ${folder_to_mount}:/data ibp-pafgrs.sif R --vanilla
```


