process query {
  tag { "${query.baseName}" }
  // The cluster's 1 GB default OOM-kills the larger queries.
  memory { 16.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 2
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(max_count)
  path(query)

  output:
  path("${query.baseName}.parquet")

  script:
  // max_count is no longer used inside the script (there's no dense-fill
  // padding step anymore - metadata-build's left join handles missing ids
  // natively) but stays as an input to preserve the dependency edge on
  // urs_counts at the call site (basic_query(urs_counts, basic_sql), etc.)
  // without having to touch every call.
  """
  rnac precompute extract-query $query ${query.baseName}.parquet
  """
}
