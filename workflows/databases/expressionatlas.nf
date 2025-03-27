process fetch_taxids {
  queue 'datamover'
  memory '8GB'

  input:
    tuple val(flag), path(experiments_path)

  output:
    path('taxids_to_fetch')

  """
  rnac expressionatlas get-taxids ${experiments_path} taxids_to_fetch
  """
}

process find_experiments {
  queue 'datamover'

  input:
    tuple val(flag), path(experiments_path)

  output:
    path('experiment_list')

  """
  find `readlink ${experiments_path}`/ -maxdepth 1 -name 'E*' -type d ! -empty > experiment_list
  """
}

process synchronize_cache {
  queue 'datamover'
  errorStrategy { task.exitStatus == 23 ? 'ignore' : 'terminate' }


  input:
    tuple path(experiments_path), path(ea_cache_path)
  output:
    val('cache synchronized')

  script:
  """
  find experiments/ -maxdepth 1 -type d ! -readable -o -type d ! -executable | sort -u >> exclude_dirs

  rsync -qLrtvz --ignore-errors --include "*analytics.tsv" \
  --filter="+ */" \
  --filter="+ *condensed-sdrf.tsv" \
  --filter="+ *-tpms.tsv" \
  --filter="+ *-configuration.xml" \
  --filter="- *-transcripts-tpms.tsv" \
  --filter="- *" \
  --exclude-from=exclude_dirs \
  ${experiments_path}/ ${ea_cache_path}/ || true
  """
}

process fetch_lookup {
  memory '8GB'
  input:
    tuple path(query), path(taxids_to_fetch)

  output:
    path("lookup_dump.csv")

    """
    psql -f $query $PGDATABASE > lookup_dump.csv
    """
}


process parse_tsvs {
  memory '2GB'
  cpus 4

  input:
  tuple path(tsvs), path(lookup)

  output:
  path('chunk_*')

  """
  rnac expressionatlas parse-dir ${tsvs} ${lookup} genes_hit.ndjson
  """

}

process group_and_convert {
  memory { 32.GB * task.attempt }
  cpus 16
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  input:
    tuple path(lookup), path(genes)

  output:
    path('*.csv')

  """
  rnac expressionatlas parse ${genes} ${lookup} .

  """
}


workflow expressionatlas {

  emit: data
  main:

  if( params.databases.expressionatlas.run ) {
    Channel.fromPath('files/import-data/expressionatlas/lookup-dump-query.sql') | set { lookup_sql }
    Channel.of(params.databases.expressionatlas.remote) | set { tsv_path }
    Channel.of(params.databases.expressionatlas.cache) | set { ea_cache }

    tsv_path.combine(ea_cache) | synchronize_cache |  set { cache_syncd }

    cache_syncd.combine(ea_cache) | fetch_taxids | set { taxids }
    cache_syncd.combine(ea_cache) | find_experiments | set { experiments }

    lookup_sql.combine(taxids) | fetch_lookup | set { lookup }

    experiments \
    | splitText \
    | combine(lookup) \
    | parse_tsvs \
    | collectFile \
    | group_and_convert \
    | collect \
    | set {data}

  }
  else {
    Channel.empty() | set { data }
  }

}


workflow {
  expressionatlas()
}
