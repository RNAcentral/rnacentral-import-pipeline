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


# store PERL5LIB
perl5lib_backup=$PERL5LIB
# read config, change PERL5LIB
. config/hive_params
# delete old output data
rm -Rf $DATA_OUT/*
rm -f log/rnacentral_import.log;
# launch the script
perl scripts/rnac_loader.pl -in=$DATA_IN \
                            -out=$DATA_OUT \
                            -user=$ORACLE_USER \
                            -password=$ORACLE_PASSWORD \
                            -host=$ORACLE_HOST \
                            -port=$ORACLE_PORT \
                            -sid=$ORACLE_SID \
                            -release_type=$DB_RELEASE_TYPE;
# restore PERL5LIB
export PERL5LIB=$perl5lib_backup
