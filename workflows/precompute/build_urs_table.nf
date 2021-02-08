include { fetch_release_info as precompute_releases } from './utils'
include { fetch_release_info as xref_releases } from './utils'

process create_schema {
  executor 'local'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(sql)

  output:
  val('done')

  """
  psql -v ON_ERROR_STOP=1 -f $sql $PGDATABASE
  """
}

process fetch_all_urs_taxid {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(_flag)
  path(query)

  output:
  path('data.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > data.csv
  """
}

process select_outdated {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory '10GB'

  input:
  path('xref.csv')
  path('precompute.csv')
  path('active.txt')

  output:
  path('to-load.csv')

  """
  precompute select xref.csv precompute.csv urs.csv
  expand-urs text active.txt urs.csv to-load.csv
  """
}

process run_query {
  tag { "$query.baseName" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(_flag)
  path(query)

  output:
  path('ids.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > ids.csv
  """
}

process build_table {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple path('computed*.csv'), path(load)

  output:
  val('done')

  """
  sort -u computed*.csv > to-load.csv
  psql \
    -v ON_ERROR_STOP=1 \
    -f "$load" "$PGDATABASE"
  """
}

workflow using_release {
    take: flag
    emit: selected
    main:
      Channel.fromPath('files/all-active-urs-taxid.sql') | set { active_sql }
      Channel.fromPath('files/precompute/fetch-xref-info.sql') | set { xref_sql }
      Channel.fromPath('files/precompute/fetch-precompute-info.sql') | set { pre_sql }

      precompute_releases(flag, pre_sql) | set { precompute_info }
      xref_releases(flag, xref_sql) | set { xref_info }
      fetch_all_urs_taxid(flag, active_sql) | set { active_urs }

      select_outdated(xref_info, precompute_info, active_urs) | set { selected }
}

workflow using_query {
  take: flag
  emit: selected
  main:
    Channel.fromPath("$params.precompute.select.query") | set { to_select }

    run_query(flag, to_select) | set { selected }
}

workflow build_urs_table {
    take: method
    emit: finished
    main:
      Channel.fromPath('files/precompute/schema.sql') | set { schema_sql }
      Channel.fromPath('files/precompute/load-urs.sql') | set { load_sql }

      create_schema(schema_sql) \
      | combine(method) \
      | map { _flag, method -> method } \
      | branch {
        release: it == 'release'
        query: it == 'query'
      } | set { to_build }

      to_build.release | using_release | set { from_release }
      to_build.query | using_query | set { from_query }

      from_release.mix(from_query, from_release) \
      | collect \
      | combine(load_sql) \
      | build_table \
      | set { finished }
}
