# Install all packages, functions, and procedures using the automatically generated
# script install_all_code.sql from the sql folder

# find directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
STARTING_DIR=`pwd`
cd $DIR
cd ../sql

# load parameters
. ../config/hive_params
echo "@install_all_code.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(HOST='$ORACLE_HOST')(PORT='$ORACLE_PORT'))(CONNECT_DATA=(SERVER=DEDICATED)(SERVICE_NAME='$ORACLE_SID')))'

# return to the original location
cd $STARTING_DIR
