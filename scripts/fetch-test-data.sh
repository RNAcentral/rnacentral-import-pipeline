#!/usr/bin/env bash

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
