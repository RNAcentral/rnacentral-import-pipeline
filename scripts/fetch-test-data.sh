#!/usr/bin/env bash

# Copyright [2009-2017] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -euo pipefail
IFS=$'\n\t'

get_data()
{
  filename="$1"
  expected="$2"
  test_file="$(basename $1)"
  final="$(basename "$filename" ".gz")"
  pushd data
  wget -O - "$filename" | gzip -d > "$final"
  hash="$(md5 -q "$final")"
  [ "$hash" = "$expected" ] || {
    echo 2>&1 "Unexpected hash for $filename";
    echo 2>&1 "Expected: $expected"
    echo 2>&1 "Actual: $hash"
    exit 1
  }
  popd
}

get_data 'ftp://ftp.ensembl.org/pub/release-87/embl/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.87.chromosome.IV.dat.gz' '2b693489f7fb5cf7db1b9f7cdaf14c50'
get_data 'ftp://ftp.ensembl.org/pub/release-87/embl/homo_sapiens/Homo_sapiens.GRCh38.87.chromosome.12.dat.gz' '51bbbf2516669d3b5fa39a2237833613'
