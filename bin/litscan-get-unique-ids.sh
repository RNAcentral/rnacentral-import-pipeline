#!/bin/bash

# set parameters
file=$1
database=$2

awk -F'|' '{
    if (NF == 2) {
        job_id = $1;
        urs = $2;
    } else {
        job_id = $1;
        primary_id = $2;
        urs = $3;
    }

    if (job_id) print job_id;
    if (primary_id) print primary_id;
    if (urs) print urs;
}' "$file" | LC_ALL=C sort -u > "${database}_ids.txt"
