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


# import tmrna Website accession data into the rnc_accessions table

# look for the file in the same directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

sqlplus -s $ORACLE_USER/$ORACLE_PASSWORD@$ORACLE_SID << EOF
whenever sqlerror exit sql.sqlcode;
set echo off
set heading off

@$DIR/update_tmrna_accessions.sql

exit;
EOF

