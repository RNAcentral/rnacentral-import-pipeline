# RNAcentral data import pipeline

## About

rnac-loader is the main RNAcentral pipeline that loads the ENA Non-coding Product into the RNAcentral Oracle staging table.

For more information, see http://www.ebi.ac.uk/seqdb/confluence/display/RNAC/RNAcentral+data+import+pipeline

## Installation

### Perl dependencies

-   BioPerl
-   DBI
-   DBD::Oracle
-   Log4perl

### eHive and PERL5LIB

	# choose the installation location
	export RNAC=/path/to/the/installation/location
	cd $RNAC

	# install ensembl hive
	mkdir $RNAC/ensembl_main
	cd $RNAC/ensembl_main

	# log in to the CVS server (password CVSUSER)
	cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login

	# checkout ensembl code
	cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout ensembl
	cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout ensembl-hive

	# set up environmental variables to enable ensembl
	export ENSEMBL_CVS_ROOT_DIR="$RNAC/ensembl_main"
	export PERL5LIB=${PERL5LIB}:${ENSEMBL_CVS_ROOT_DIR}/ensembl/modules
	export PERL5LIB=${PERL5LIB}:${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/modules

	# enable RNAcentral modules
	export PERL5LIB=${PERL5LIB}:${RNAC}/rnac-loader/modules;
	# enable RNAcentral-specific Hive modules
	export PERL5LIB=${RNAC}/rnac-loader/lib/:${PERL5LIB};

	# clone rnac-loader
	cd $RNAC
	git clone https://apetrov@scm.ebi.ac.uk/git/rnac-loader.git

	# create an empty mysql Hive database
	mysql -h <hostname> -u <user> --port <port> -p
	> CREATE DATABASE IF NOT EXISTS rnacentral_staging;

## Configuration

The options are specified as command line arguments and as default options in the Bio::RNAcentral::Base module.

	 cp config/hive_params_template config/hive_params
	# add sensitive connection details to config/hive_params
	. config/hive_params

## Running in Hive mode

	. scripts/hive_pipeline.sh

## Running in non-Hive mode

	perl scripts/rnac_loader.pl <options>

## License

	See the LICENSE file for more information.

