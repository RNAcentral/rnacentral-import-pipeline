# Install all packages, functions, and procedures into a single Oracle database
# using the automatically generated script install_all_code.sql from the sql folder

# find directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
STARTING_DIR=`pwd`
cd $DIR
cd ../sql

# load parameters
. ../config/hive_params
echo "@install_all_code.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@'$ORACLE_SID

# return to the original location
cd $STARTING_DIR
