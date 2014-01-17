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


#################
## Data import ##
#################
# store PERL5LIB
perl5lib_backup=$PERL5LIB
# read config, change PERL5LIB
. config/hive_params
# delete old output data
rm -Rf $DATA_OUT/*
rm -f log/rnacentral_import.log;
# initialize hive pipeline
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::RNAcentralUpdate_conf \
	-output_folder=$DATA_OUT \
	-release_type=$DB_RELEASE_TYPE \
	-pipeline-db -pass=$HIVE_PASSWORD \
	-pipeline-db -user=$HIVE_USERNAME \
	-pipeline-db -host=$HIVE_HOST \
	-pipeline-db -port=$HIVE_PORT \
	-oracle-user=$ORACLE_USER \
	-oracle-password=$ORACLE_PASSWORD \
	-oracle-host=$ORACLE_HOST \
	-oracle-sid=$ORACLE_SID \
	-oracle-port=$ORACLE_PORT;
HIVE_URL='mysql://'$HIVE_USERNAME':'$HIVE_PASSWORD'@'$HIVE_HOST':'$HIVE_PORT'/rnacentral_staging'
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/beekeeper.pl -url $HIVE_URL -sync;
# generate graph
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/generate_graph.pl -url $HIVE_URL -output pipeline.png
# launch the pipeline
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/beekeeper.pl -url $HIVE_URL -meadow_type LSF -loop -total_running_workers_max 40;

# manually import tmrna website accessions
. scripts/tmrna_website_manual_import/update_tmrna_accessions.sh

###############
## Reporting ##
###############
# generate title
if grep --quiet "WARN\|FATAL" log/rnacentral_import.log; then
    title='[FAILED] RNAcentral data import'
else
    title='[OK] RNAcentral data import'
fi
# generate body
echo 'http://test.rnacentral.org' > email_message_body.txt
grep "WARN\|FATAL" log/rnacentral_import.log >> email_message_body.txt
# compress log
gzip log/rnacentral_import.log
# email the report
mutt $RNACENTRAL_ADMIN_EMAIL -s "$title" -a log/rnacentral_import.log.gz -i email_message_body.txt < /dev/null
# remove body
rm email_message_body.txt
# uncompress log
gzip -d log/rnacentral_import.log.gz

# restore PERL5LIB
export PERL5LIB=$perl5lib_backup
