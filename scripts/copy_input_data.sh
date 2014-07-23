#!/bin/bash
# Copyright [2009-2014] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


###############################################################
## Copy input text files to another directory for processing ##
###############################################################

# The input files are found in several location (Non-coding, RFAM, RefSeq  products).
# This script copies the files into a location accessible from the LSF farm.
#
# Usage:
# . copy_input_data.sh /path/to/input/product /path/to/processing/location

STARTING_DIR=`pwd`

rnacentral_product=$1
final_location=$2

cd $rnacentral_product

# dat files
for files in *.dat
do
	cp "$files" "${final_location}/${files%.dat}.ncr"
done

# ncr files
for files in *.ncr
do
	cp "$files" "${final_location}/$files"
done

# ncr.gz files
for files in `find . -name '*.ncr.gz'`
do
	# $files - nested file path, $filename - file basename
	$filename=`basename $files`
	cp "$files" "${final_location}/$filename"
	gzip -d "${final_location}/$filename"
done

# return to the original location
cd $STARTING_DIR
