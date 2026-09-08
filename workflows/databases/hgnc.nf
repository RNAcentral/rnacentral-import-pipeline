process hgnc {
  output:
  path('*.{csv,parquet}')

  when: params.databases.hgnc?.run

  script:
  def force = params.force_full_import ? '--force-full' : ''
  """
  wget -O raw.json $params.databases.hgnc.remote
  rnac hgnc map $force raw.json
  """
}
