#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

name=$(git rev-parse --abbrev-ref HEAD)
sif="$name.sif"

tag="$name"
if [[ "$name" = "master" ]]; then
  tag="latest"
  sif="main.sif"
fi

# Build docker image
cp ../requirements.txt .
docker build -t "bsweeneyebi/rnacentral-import-pipeline:$tag" .
docker push "bsweeneyebi/rnacentral-import-pipeline:$tag"

# Build singuarlity image
vagrant up
vagrant ssh -c "sudo singularity build $sif docker://bsweeneyebi/rnacentral-import-pipeline:$tag"
vagrant ssh -c "cat /home/vagrant/$sif" > "$sif"
vagrant halt
chmod +x "$sif"

# scp main.sif "$USER@yoda-login:/nfs/leia/production/xfam/users/bsweeney/containers/"
