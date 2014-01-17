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


# Install all packages, functions, and procedures into a single Oracle database
# using the automatically generated script install_all_code.sql from the sql folder

# find directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
STARTING_DIR=`pwd`
cd $DIR
cd ../sql

# load parameters
. ../config/hive_params
echo "@install_all_code.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@'$ORACLE_SID

# return to the original location
cd $STARTING_DIR
