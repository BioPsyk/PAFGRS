#!/usr/bin/env bash

if [ -z "${folder_to_mount}" ]; then
    folder_to_mount="$HOME"
fi
echo "folder_to_mount is set to: $folder_to_mount, destination is /data"

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${script_dir}/init-containerization.sh"

mount_flags=$(format_mount_flags "-B")

FAKE_HOME="tmp/fake-home"
export SINGULARITY_HOME="/pafgrs/${FAKE_HOME}"
mkdir -p "${FAKE_HOME}"

exec singularity run \
     --contain \
     --cleanenv \
     ${mount_flags} \
     --bind ${folder_to_mount}:/data \
     "tmp/${singularity_image_tag}" \
     "$@"
