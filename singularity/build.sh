#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

function cleanup() {
  rm -f "$ssh_config"
}

name="$(git rev-parse --abbrev-ref HEAD)"
sif="$name.sif"

tag="$name"
if [[ "$name" = "master" ]]; then
  tag="latest"
fi

docker="rnacentral/rnacentral-import-pipeline:$tag"
ssh_config="$(mktemp)"

cd "$(dirname $0)"

vagrant up
vagrant ssh-config > "$ssh_config"
echo ssh -F "$ssh_config" default "sudo singularity build $sif docker://$docker"
ssh -F "$ssh_config" default "sudo singularity build $sif docker://$docker"

# ssh -F "$ssh_config" default "su -; sudo singularity build $sif docker://$docker"
# scp -F "$ssh_config" -P 2222 "vagrant@127.0.0.1:/home/vagrant/$sif" .
# vagrant halt

trap cleanup EXIT
