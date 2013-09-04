. config/hive_params
rm -Rf $DATA_OUT/*
rm -f log/rnacentral_import.log;
perl scripts/rnac_loader.pl -in=$DATA_IN \
                            -out=$DATA_OUT \
                            -user=$ORACLE_USER \
                            -password=$ORACLE_PASSWORD \
                            -host=$ORACLE_HOST \
                            -port=$ORACLE_PORT \
                            -sid=$ORACLE_SID \
                            -release_type=$DB_RELEASE_TYPE;
