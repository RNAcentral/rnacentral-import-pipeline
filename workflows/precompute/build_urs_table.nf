include { fetch_release_info as precompute_releases } from './utils'
include { fetch_release_info as xref_releases } from './utils'

process fetch_all_urs_taxid {
  when { params.precompute.run }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(query)

  output:
  path('data.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > data.csv
  """
}

process build_table {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory '10GB'

  input:
  tuple path(load), path('xref.csv'), path('precompute.csv'), path('active.txt')

  output:
  val('done')

  """
  precompute select xref.csv precompute.csv urs.csv
  expand-urs text active.txt urs.csv to-load.csv
  psql \
    --variable ON_ERROR_STOP=1 \
    --variable tablename=${params.precompute.tablename} \
    -f "$load" "$PGDATABASE"
  """
}

workflow build_urs_table {
    emit: finished
    main:
      Channel.fromPath('files/precompute/load-urs.sql') | set { load_sql }
      Channel.fromPath('files/all-active-urs-taxid.sql') | set { active_sql }
      Channel.fromPath('files/precompute/fetch-xref-info.sql') | set { xref_sql }
      Channel.fromPath('files/precompute/fetch-precompute-info.sql') | set { pre_sql }

      precompute_releases(pre_sql) | set { precompute_info }
      xref_releases(xref_sql) | set { xref_info }
      fetch_all_urs_taxid(active_sql) | set { active_urs }

      load_sql \
      | combine(xref_info) \
      | combine(precompute_info) \
      | combine(active_urs) \
      | build_table \
      | set { finished }
}
