# How to run an RNAcentral release

This outlines all the steps needed to run a release. This doesn't cover any of
the debugging that may be needed if a database changes how they provide data.

## Setup environment

There are some required directories and environment variables that have to be
setup for the pipeline to run.

### Clone the code

First clone the code into some directory. I do:

```sh
cd /hps/nobackup/production/xfam/bsweeney
git clone RNAcentral/rnacentral-import-pipeline.git release-number
cd release-number
git checkout dev
```

### Compile the code

The pipeline uses some rust code which needs to be compiled ahead of time. Do
this with:

```sh
$ make rust
```

This puts the programs into `bin/` which will be part of `$PATH` inside
nextflow.

### Create required directories

We use a temp directory inside of `work` to try and avoid a possibly slow
filesystem on nodes. 

```
mkdir -p work/tmp
```

### Set required environment variables

Set the required env variables.

```sh
$ export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'
```

This variable is needed to ensure nextflow can create enough threads for the
large number of jobs the pipeline runs.

### Download required R2DT files

Fetch the latest R2DT files from the ftp site and place them where they are
expected. Note, you may need to update the version `1.2` below to the latest
one.

```sh
$ mkdir -p singularity/bind/r2dt/data/cms
$ cd singularity/bind/r2dt/data/cms
$ wget 'http://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/1.2/cms.tar.gz'
$ tar xvf cms.tar.gz
```

### Optionally install nextflow

You may also need to install nextflow depending on what cluster you are using.
This can be done with:

```sh
$ curl -s https://get.nextflow.io | bash
```

which will create a `./nextflow` script. In all follow commands use
`./nextflow` instead of `nextflow`.

## Configure the pipeline

The pipeline is configured through a file called `local.config`. This file
controls what databases are imported and what steps are run, eg will Rfam
analysis run or not. This is done because it is common to need to rerun some
steps and by setting a simple toggle this will turn off attempting to run
certain steps.

Here is an example:

```groovy
params {
  notify = false
  release = 19

  databases {
    ena.run = true
    ensembl {
      plants.run = true
      fungi.run = true
      protists.run = true
      metazoa.run = true
      vertebrates.run = true
    }
    flybase.run = true
    genecards.run = true
    gtrnadb.run = true
    hgnc.run = true
    intact.run = true
    lncipedia.run = true
    malacards.run = true
    pdbe.run = true
    pombase.run = true
    quickgo.run = true
    refseq.run = true
    sgd.run = true
    silva.run = true
    zfin.run = true
  }

  genome_mapping {
    run = true
    select_mapped.directives.memory = 8.GB
    blat.directives.memory = 10.GB
    downloads.directives.memory = 12.GB
  }

  cpat.run = true

  qa {
    rfam.run = false
    rfam.memory = 8.GB
    dfam.run = false
    pfam.run = false
  }

  precompute {
    run = true
    maxForks = 4
    range.memory = '5GB'
    method = 'query'
    select.query = 'changed-models.sql'
  }

  r2dt {
    run = false
    publish = "/nfs/production/xfam/rnacentral/secondary-structure"
  }

  feedback.run = false

  search_export {
    max_entries = 50000
    max_forks = 1
    memory = '15 GB'
    publish {
      host = ''
      path = "/nfs/production/xfam/rnacentral/search_dumps/dev-nightly/"
    }
  }

  sequence_search {
    run = true
  }
}

singularity {
  enabled = true
  cacheDir = "$baseDir/singularity"
  runOptions = '--bind /nfs/ftp/pub/databases/ena --bind /ebi/ftp --bind /nfs/ftp --bind /nfs/ensemblftp --bind /nfs/ensemblgenomes/ftp'
}

includeConfig '../profiles.config'
includeConfig 'config/cluster.config'
```

## Import data

Run the workflow to import data with:

```sh
$ nextflow run -profile prod import-data.nf
```

## Analyze sequences

This steps does:

- Rfam scan
- R2DT scan
- CPAT scan
- Genome mapping

These can all be done with:

```sh
$ nextflow run -profile prod analyze.nf
```

## Precompute

This computes the RNA types, descriptions and QA analysis. It is run with:

```sh
$ nextflow run -profile prod precompute.nf
```

## Build genes

This builds RNA genes and must come after precompute as it depends on the RNA
types of individual sequences.

```sh
$ nextflow run -profile prod genes.nf
```

## Export data

This builds the text search, sequence search and ftp exports. It is run with:

```sh
$ nextflow run -profile prod export.nf
```

## Other notes

There is a script `v2.nf` that is meant to run all the workflows in order. It
is untested currently.
