# RNAcentral data import pipeline

## About

This is the main pipeline that is used internally for loading the data into the
RNAcentral database. [More
information](http://www.ebi.ac.uk/seqdb/confluence/display/RNAC/RNAcentral+data+import+pipeline).
The pipeline is [nextflow](https://www.nextflow.io) based and the main entry
point is main.nf. 

The pipeline is typically run as:

```sh
nextflow run -profile env -with-singularity pipeline.sif main.nf
```

The pipeline is meant to run 

## Configuring the pipeline

The pipeline requires a `local.config` file to exist and contain some
information. Notably a `PGDATABASE` environment variable must be defined so
data can be imported or fetched. In addition, to import specific databases
there must be a `params.import_data.databases` dict defined. The keys must be
known databases names and the values should be truthy to indicate the databases
should be imported.

There is some more advanced configuration options available, such as turning on or off
specific parts of the pipeline like genome mapping, qa, etc.

## Using with Docker

The pipeline is meant to run in docker or singularity. You should build or
fetch a suitable container. Some example commands are below.

* build container
  ```
  docker build -t rnacentral-import-pipeline .
  ```

* open interactive shell inside a running container
  ```
  docker run -v `pwd`:/rnacentral/rnacentral-import-pipeline -v /path/to/data:/rnacentral/data/ -it rnacentral-import-pipeline bash
  ```

## Testing

Several tests require fetching some data files prior to testing. The files can
be fetched with:

```sh
./scripts/fetch-test-data.sh
```

The tests can then be run using [py.test](http://pytest.org). For example,
running Ensembl importing tests can be done with:

```sh
py.test tests/databases/ensembl/
```

## Other environment variables

The pipeline requires the `NXF_OPTS` environment variable to be set to
`-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000`, a module for doing this is
in `modules/cluster`. Also some configuration settings for efficient usage on
EBI's LSF cluster  are in `config/cluster.config`.

## License

See [LICENSE](https://github.com/RNAcentral/rnacentral-import-pipeline/blob/master/LICENSE) for more information.
