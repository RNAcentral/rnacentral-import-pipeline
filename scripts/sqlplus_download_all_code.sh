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


# Download all procedures, packages, and functions from the Oracle database
# and store them as individual files in the same directory as this script

# find directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
STARTING_DIR=`pwd`
cd $DIR
cd ../sql

# load parameters
. ../config/hive_params
echo "@getallcode.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(HOST='$ORACLE_HOST')(PORT='$ORACLE_PORT'))(CONNECT_DATA=(SERVER=DEDICATED)(SERVICE_NAME='$ORACLE_SID')))'
rm xtmpx.sql

# return to the original location
cd $STARTING_DIR
