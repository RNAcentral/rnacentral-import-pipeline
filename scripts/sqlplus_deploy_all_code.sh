# Deploy all packages, functions, and procedures to all Oracle instances
# using the automatically generated script install_all_code.sql from the sql folder.

# find directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
STARTING_DIR=`pwd`
cd $DIR
cd ../sql

echo "@install_all_code.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@EGRNAPRO'
echo "@install_all_code.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@EGRNADEV'
echo "@install_all_code.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@EGRNATST'
echo "@install_all_code.sql" | sqlplus $ORACLE_USER'/'$ORACLE_PASSWORD'@EGRNARLS'

# return to the original location
cd $STARTING_DIR
