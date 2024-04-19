#!/usr/bin/env bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
project_dir=$(dirname "${script_dir}")

function format_mount_flags() {
  flag="${1}"

  for mount in "${mounts[@]}"
  do
    echo "${flag} ${project_dir}/${mount}:/pafgrs/${mount} "
  done
}

cd "${project_dir}"

# mount the whole project root
mounts=(
  $(ls -a -1 $PWD | tail -n+3)
)

image_tag="ibp-pafgrs-base:"$(cat "docker/VERSION")
deploy_image_tag="ibp-pafgrs:"$(cat "docker/VERSION")

#singularity build
singularity_image_tag="ibp-pafgrs-base_version-$(cat "docker/VERSION").sif"

mkdir -p tmp

