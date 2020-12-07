#!/bin/bash

# set the file based on the first argument and the environment based on the second argument
file=$1
env=$2

# set the name of the file to be copied and the S3 path
urs=$(echo "$file" | grep -o 'URS0[^[:blank:]]*')
path=$env/${urs:0:3}/${urs:3:2}/${urs:5:2}/${urs:7:2}/${urs:9:2}/

# name of the S3 bucket
bucket_name="ebi-rnacentral"

# get the current date to calculate the signature and also to pass to S3
date=`date +'%a, %d %b %Y %H:%M:%S %z'`

# calculate the signature to be sent as a header
content_type="application/octet-stream"
string_to_sign="PUT\n\n$content_type\n${date}\n/${bucket_name}/${path}${urs}"
signature=$(echo -en "${string_to_sign}" | openssl sha1 -hmac "${SECRET}" -binary | base64)

# upload file
echo "Adding ${urs} to S3"
curl -X PUT -T "${file}" \
     -H "Host: s3.embassy.ebi.ac.uk/${bucket_name}" \
     -H "Date: $date" \
     -H "Content-Type: $content_type" \
     -H "Authorization: AWS ${S3_KEY}:${signature}" \
     "https://s3.embassy.ebi.ac.uk/${bucket_name}/${path}${urs}"
