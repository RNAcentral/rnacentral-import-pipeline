#!/bin/bash
. config/hive_params
rm log/rnacentral_import.log;
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
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/generate_graph.pl -url $HIVE_URL -output pipeline.png
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/beekeeper.pl -url $HIVE_URL -meadow_type LSF -loop -total_running_workers_max 40;
