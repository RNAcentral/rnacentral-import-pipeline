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


########################################################
## Add RFAM product files to Non-coding product files ##
########################################################

# RFAM product resides in a separate location from the Non-coding product
# This script moves the RFAM files into the same location
# as the rest of the Non-coding product.
# The idea is to run this script before a full update or when the RFAM product
# needs to be updated.
#
# Usage:
# . copy_rfam_product.sh /path/to/rfam/product /path/to/noncoding/product

STARTING_DIR=`pwd`

rfam_product=$1
noncoding_product=$2

cd $rfam_product

for files in *.dat
do
	cp "$files" "${noncoding_product}/${files%.dat}.ncr"
done

# return to the original location
cd $STARTING_DIR
