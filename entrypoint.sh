#!/bin/sh
source $LOC/virtualenvs/rnacentral/bin/activate
cd /rnacentral/rnacentral-import-pipeline
exec "$@"
