# RNAcentral data import pipeline

## About

This is the main pipeline that is used internally for loading the data into the RNAcentral database.
[More information](http://www.ebi.ac.uk/seqdb/confluence/display/RNAC/RNAcentral+data+import+pipeline)

## Installation

### Perl dependencies

-   BioPerl
-   DBI
-   DBD::Oracle
-   Log4perl

### Installation

```
git clone https://github.com/RNAcentral/rnacentral-import-pipeline.git
cd rnacentral-import-pipeline
# initialize submodules
git submodule init
git submodule update
# add sensitive connection details to config/hive_params
cp config/hive_params_template config/hive_params
```

## Running in Hive mode

	source scripts/hive_pipeline.sh

## Running in non-Hive mode


	perl scripts/rnac_loader.pl <options>

or

	source scripts/rnac_loader_wrapper.sh

## Testing

The tests can be run using [py.test](http://pytest.org). For example, running
Ensembl importing tests can be done with:

```sh
py.test luigi/tests/ensembl_test.py
```

## License

See [LICENSE](https://github.com/RNAcentral/rnacentral-import-pipeline/blob/master/LICENSE) for more information.
