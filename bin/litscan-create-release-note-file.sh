#!/bin/bash

# set parameter
directory=$1
version=$2

read total <<< $(zgrep '<entry id=' ${directory}/*.xml.gz | wc -l)
release_date=$(date +%d-%b-%Y)
cat > release_note.txt <<EOF
release=$version
release_date=$release_date
entries=$total
EOF
