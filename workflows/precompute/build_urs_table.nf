process create_schema {
  executor 'local'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(sql)

  output:
  val('done')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $sql \$PGDATABASE
  """
}

process select_outdated {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory '24 GB'
  cpus 4

  input:
  path('xref.csv')
  path('precompute.csv')

  output:
  path('urs.csv')

  script:
  """
  LC_ALL=C sort -t, -k1,1 --parallel=${task.cpus} -S 4G xref.csv > xref.sorted.csv
  LC_ALL=C sort -t, -k1,1 --parallel=${task.cpus} -S 4G precompute.csv > precompute.sorted.csv
  precompute select xref.sorted.csv precompute.sorted.csv urs.csv
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

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $query \$PGDATABASE > ids.csv
  """
}

process build_table {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory '10GB'

  input:
  tuple path('computed*.csv'), path(load), path('counts.sql')

  output:
  path('counts.txt')

  script:
  """
  sort -u computed*.csv > to-load-urs-taxid.csv
  test -s to-load-urs-taxid.csv || { echo "No urs_taxid selected to precompute" >&2; exit 1; }
  psql \
    -v ON_ERROR_STOP=1 \
    -f "$load" "\$PGDATABASE"
  psql -f counts.sql -v ON_ERROR_STOP=1 "\$PGDATABASE" > counts.txt
  """
}

process sort_ids {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(_flag)
  path('raw*')

  output:
  path('finalized')

  script:
  """
  sort -u raw* > finalized
  """
}

process xref_releases {
  input:
  val(_flag)
  file(query)

  output:
  path('data.csv')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $query \$PGDATABASE > data.csv
  """
}

process precompute_releases {
  input:
  val(_flag)
  file(query)

  output:
  path('data.csv')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $query \$PGDATABASE > data.csv
  """
}

workflow using_release {
    take: flag
    main:
      channel.fromPath('files/precompute/fetch-xref-info.sql') | set { xref_sql }
      channel.fromPath('files/precompute/fetch-precompute-info.sql') | set { pre_sql }

      precompute_releases(flag, pre_sql) | set { precompute_info }
      xref_releases(flag, xref_sql) | set { xref_info }

      select_outdated(xref_info, precompute_info) | set { selected }
    emit: selected
}

workflow using_query {
  take: flag
  main:
    flag \
    | map { _flag -> file(params.precompute.select.query) } \
    | set { to_select }

    run_query(flag, to_select) | set { selected }
  emit: selected
}

workflow using_all {
  take: flag
  main:
    channel.fromPath("files/precompute/methods/all.sql") | set { to_select }

    run_query(flag, to_select) | set { selected }
  emit: selected
}

workflow using_ids {
  take: flag
  main:

    flag \
    | map { _flag -> file(params.precompute.select.id_file) } \
    | set { id_files }

    sort_ids(flag, id_files) | set { selected }
  emit: selected
}

workflow build_urs_table {
    take: method
    main:
      channel.fromPath('files/precompute/schema.sql') | set { schema_sql }
      channel.fromPath('files/precompute/load-urs.sql') | set { load_sql }
      channel.fromPath('files/precompute/get-urs-count.sql') | set { count_sql }

      create_schema(schema_sql) \
      | combine(method) \
      | map { _flag, method_val -> method_val } \
      | branch { v ->
        release: v == 'release'
        query: v == 'query'
        all: v == 'all'
        ids: v == 'ids'
      } \
      | set { to_build }

      to_build.release | using_release | set { from_release }
      to_build.query | using_query | set { from_query }
      to_build.all | using_all | set { from_all }
      to_build.ids | using_ids | set { from_ids }

      channel.empty() \
      | mix(from_release, from_query, from_all, from_ids) \
      | collect \
      | combine(load_sql) \
      | combine(count_sql) \
      | build_table \
      | splitCsv \
      | first \
      | map { row -> row[0].toInteger() + 1 } \
      | set { finished }
    emit: finished
}
