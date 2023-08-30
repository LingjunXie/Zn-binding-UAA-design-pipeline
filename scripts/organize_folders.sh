#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -uaa) uaa_name="$2";;
        *) break;
    esac
    shift 2
done

# Create a new folder to move the folders into
mkdir ${uaa_name}

# Loop through all folders and move those ending with "APC" to the new folder
for folder in *_${uaa_name}; do
  if [ -d "$folder" ]; then
    mv "$folder" ${uaa_name}/
  fi
done
