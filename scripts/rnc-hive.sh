perl ../ensembl-hive/scripts/init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::RNAcentralUpdate_conf -pipeline-db -host=$HOST -pipeline-db -user=$USER -pipeline-db -pass=$PASSWORD;
perl ../ensembl-hive/scripts/beekeeper.pl -url mysql://$USER:$PASSWORD@$HOST:3306/apetrov_rnacentral_update -sync;
perl ../ensembl-hive/scripts/beekeeper.pl -url mysql://$USER:$PASSWORD@$HOST:3306/apetrov_rnacentral_update -meadow_type LOCAL -total_running_workers_max 5 -loop;
