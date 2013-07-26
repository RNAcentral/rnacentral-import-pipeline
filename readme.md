# RNAcentral data import pipeline

## About

rnac-loader is the main RNAcentral pipeline that loads the ENA Non-coding Product into the RNAcentral Oracle staging table.

## Installation

### Perl dependencies

-   BioPerl
-   DBI
-   DBD::Oracle
-   Log4perl

### eHive

	$ mkdir $HOME/ensembl_main

	$ cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login
	Logging in to :pserver:cvsuser@cvs.sanger.ac.uk:2401/cvsroot/ensembl
	CVS password: CVSUSER

	$ cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout ensembl
	$ cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout ensembl-hive

	$ export ENSEMBL_CVS_ROOT_DIR="$HOME/ensembl_main"
	$ export PERL5LIB=${PERL5LIB}:${ENSEMBL_CVS_ROOT_DIR}/ensembl/modules
	$ export PERL5LIB=${PERL5LIB}:${ENSEMBL_CVS_ROOT_DIR}/ensembl-hive/modules

	# enable RNAcentral modules
	$ export PERL5LIB=${PERL5LIB}:${HOME}/rnac-loader/modules;
	# enable RNAcentral-specific Hive modules
	$ export PERL5LIB=${HOME}/rnac-loader/lib/:${PERL5LIB};

## Configuration

The options are specified as command line arguments and as default options in the Bio::RNAcentral::Base module.

## Running in Hive mode

	. ./scripts/rnac-hive.sh

## Running in non-Hive mode

	perl scripts/rnac_loader.pl <options>

