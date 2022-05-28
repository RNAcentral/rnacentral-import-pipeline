#!/bin/bash

# set parameters
file=$1
database=$2

# read file line by line
while IFS= read -r line; do
    IFS=$"|"
    tmp=($line)
    if [[ ${#tmp[*]} = 2 ]]; then
      job_id="${tmp[0]}"
      urs="${tmp[1]}"
    else
      job_id="${tmp[0]}"
      primary_id="${tmp[1]}"
      urs="${tmp[2]}"
    fi

    if [[ -n "${job_id}" ]]; then
      echo ${job_id} >> ${database}_all_ids.txt
    fi

    if [[ -n "${primary_id}" ]]; then
      echo ${primary_id} >> ${database}_all_ids.txt
    fi

    if [[ -n "${urs}" ]]; then
      echo ${urs} >> ${database}_all_ids.txt
    fi
done < ${file}

# create file with unique ids
cat ${database}_all_ids.txt | sort | uniq > ${database}_ids.txt
