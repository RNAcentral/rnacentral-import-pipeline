#!/bin/bash

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
