#!/usr/bin/env bash

if [ -z "${folder_to_mount}" ]; then
    folder_to_mount="$HOME"
fi
echo "folder_to_mount is set to: $folder_to_mount, destination is /data"

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${script_dir}/init-containerization.sh"

mount_flags=$(format_mount_flags "-v")

exec docker run -it --rm ${mount_flags} -v ${folder_to_mount}:/data "${image_tag}" "$@"
