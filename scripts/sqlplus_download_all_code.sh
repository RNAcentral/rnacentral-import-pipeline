# Download all procedures, packages, and functions from the Oracle database
# and store them as individual files in the same directory as this script

# find directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
STARTING_DIR=`pwd`
cd $DIR
cd ../sql

# load parameters
. ../config/hive_params
echo "@getallcode.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(HOST='$ORACLE_HOST')(PORT='$ORACLE_PORT'))(CONNECT_DATA=(SERVER=DEDICATED)(SERVICE_NAME='$ORACLE_SID')))'
rm xtmpx.sql

# return to the original location
cd $STARTING_DIR
