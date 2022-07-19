#!/bin/bash

# Script to submit ids to RNAcentral-reference

# set parameter
file=$1

# create folder
[ ! -d submitted ] && mkdir submitted

function submitJob
{
  # set parameter
  job_id=$1

  # submit job
  curl -X POST \
       -H "Content-Type:application/json" \
       -d "{\"id\": \"${job_id}\"}" \
       http://45.88.80.122:8080/api/submit-job && echo ${job_id} >> submitted/${file};
}

# loop through the file
while IFS="" read -r line || [ -n "$line" ]
do
  submitJob "$line"
done < "$file"
