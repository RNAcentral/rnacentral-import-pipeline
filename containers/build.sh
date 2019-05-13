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

docker="bsweeneyebi/rnacentral-import-pipeline:$tag"

# Build docker image
cp ../requirements.txt .
docker build -t "$docker" .
docker push "$docker"

# Build singuarlity image
vagrant up
vagrant ssh -c "su -; sudo singularity build $sif docker://$docker"

expect <<EOF
spawn scp -P 2222 vagrant@127.0.0.1:/home/vagrant/$sif . 
expect {
  password: {send "vagrant\r"; exp_continue}
}
EOF
vagrant halt

scp $sif "$USER@yoda-login:/hps/nobackup2/xfam/rnacentral/containers/"
