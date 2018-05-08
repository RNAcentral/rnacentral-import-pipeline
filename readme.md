# RNAcentral data import pipeline

## About

This is the main pipeline that is used internally for loading the data into the RNAcentral database.
[More information](http://www.ebi.ac.uk/seqdb/confluence/display/RNAC/RNAcentral+data+import+pipeline)

## Installation

## Using with Docker

* build container
  ```
  docker build -t rnacentral-import-pipeline .
  ```

* open interactive shell inside a running container
  ```
  docker run -v `pwd`:/rnacentral/rnacentral-import-pipeline -v /path/to/data:/rnacentral/data/ -it rnacentral-import-pipeline bash
  ```

* example luigi command

  ```
  python -m luigi --module tasks.release LoadRelease --local-scheduler
  ```

## Running luigi tasks

Several databases are imported using the
[luigi](https://github.com/spotify/luigi) pipeline. The code for the pipeline
is stored in `luigi` directory. The rfam search task are stored in the `tasks`
subdirectory. These can be run with:

```sh
export PYTHONPATH=$PYTHONPATH:luigi
python -m luigi --module tasks <task>
```

Sadly, luigi doesn't seem to provide a nice way to inspect the available tasks,
so the easiest way to see what is available is to read
`luigi/tasks/__init__.py`.

Some individual examples are:

```sh
python -m luigi --module tasks RfamCSV
python -m luigi --module tasks RfamSearches
```

For details on each individual part read the documentation for the task you are
interested in.

There are also several other database, like NONCODE and Greengenes, that aren't
yet moved into the tasks directory. These can be found under the `luigi/`
directory. Running these is similar, some examples are:

```sh
python -m luigi --module json_batch_processor Noncode [options]
python -m luigi --module ensembl.species SpeciesImporter [options]
```

The pipeline requires the: `luigi.cfg` file be filled out, an example file,
with comments is in `luigi.cfg.txt`. In addition there is documentation about
the configuration in `luigi/tasks/config.py`.

## Testing

Running tests for ensembl import requires downloading data from Ensembl first.
This can be done with:

```sh
./scripts/fetch-test-data.sh
```

The tests can then be run using [py.test](http://pytest.org). For example,
running Ensembl importing tests can be done with:

```sh
py.test luigi/tests/ensembl_test.py
```

## License

See [LICENSE](https://github.com/RNAcentral/rnacentral-import-pipeline/blob/master/LICENSE) for more information.
