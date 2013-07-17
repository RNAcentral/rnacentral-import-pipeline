mysql -u $USERNAME -p$PASSWORD -e 'drop database apetrov_rnacentral_update';
perl /Users/apetrov/Desktop/ensembl_main/ensembl-hive/scripts/init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::RNAcentralUpdate_conf -pipeline-db -host=$HOST -pipeline-db -user=$USER -pipeline-db -pass=$PASSWORD ;
perl /Users/apetrov/Desktop/ensembl_main/ensembl-hive/scripts/beekeeper.pl -url mysql://$USER:$PASSWORD@$HOST:3306/apetrov_rnacentral_update -sync;
perl /Users/apetrov/Desktop/ensembl_main/ensembl-hive/scripts/beekeeper.pl -url mysql://$USER:$PASSWORD@$HOST:3306/apetrov_rnacentral_update -meadow_type LOCAL -total_running_workers_max 4 -loop;
