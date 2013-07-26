perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::RNAcentralUpdate_conf \
	-dir='/Users/apetrov/Desktop/ncrna_sample_files' \
	-out='/Users/apetrov/Desktop/ensembl_main/rnac-loader/temp' \
	-hive_driver sqlite \
	-pipeline-db -pass=$HIVE_PASSWORD \
	-oracle-user=$ORACLE_USER \
	-oracle-password=$ORACLE_PASSWORD \
	-oracle-host=$ORACLE_HOST \
	-oracle-sid=$ORACLE_SID \
	-oracle-port=$ORACLE_PORT;
HIVE_URL='sqlite:///'$HIVE_USERNAME'_rnacentral_update'
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/beekeeper.pl -url $HIVE_URL -sync;
perl $ENSEMBL_CVS_ROOT_DIR/ensembl-hive/scripts/beekeeper.pl -url $HIVE_URL -meadow_type LOCAL -total_running_workers_max 4 -loop;
