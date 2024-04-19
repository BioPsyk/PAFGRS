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
load_all()
build_vignettes()
```

The vignette in the format .html, .R and .RMD can now be found in the `doc` folder. 

## Containers
For the ability to run PAFGRS in places where new installations are not possible we provide a container image using docker and singularity. First use Docker to build a docker image on a machine where you have root access. Then use singularit to convert the docker image to a singularity image, which is a file that easily can be moved and copied. 

### Singularity at HPC:s
Sometimes it is already installed and you do not need to install it, and sometimes you are only allowed to run singularities in interactive jobs or batch jobs. And sometimes you need to load a module, see example below:

```
# To load singularity on HPC:s with module systems
module load tools
module load singularity/4.0.0
```

### Download and enter the container
The container can be downloaded from dockerhub, and ad-hoc converted to a singularity .sif image. Then you just have to enter the image, load the correct libraries. Remember to mount a parent folder that contains all input data needed for PAFRGS.
```
# Get singularity from dockerhub
mkdir -p tmp
singularity pull tmp/ibp-pafgrs-base_version-1.0.0.sif docker://biopsyk/ibp-pafgrs:1.0.0

# set folder to mount, e.g., your home folder
folder_to_mount=${HOME}

# Run the image starting R directly, mounted data accessible in /data inside the image.
./scripts/singularity-run.sh

```
PAFGRS is now loaded and ready in the interactive R session

### Re-build the container
When changing the Dockerfile, the docker and singularity images need to be rebuild. This is how you can do it locally. Requires docker and singularity to be installed.
```
# Build docker image
./scripts/docker-build.sh
# Convert to Singularity image
./scripts/singularity-build.sh
```


