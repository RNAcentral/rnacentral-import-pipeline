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


# find all ncr files that might be truncated
# requires full path as a parameter
for file in "$1"/*.ncr
do
	if tail -2 "$file" | grep -q ^\/\/
	then
	  :
	else
	  echo "No match with the pattern in $file"
	fi
done
