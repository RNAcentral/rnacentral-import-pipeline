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

