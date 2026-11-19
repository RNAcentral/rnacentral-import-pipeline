process parse_ena_coding {
  tag { source.name }
  memory '4GB'

  input:
  val source

  output:
  path('*.{csv,parquet}')

  when: params.databases.ena_coding?.run

  script:
  def assembly = source.assembly ? "--assembly '${source.assembly}'" : ''
  def end = source.end ? "--end ${source.end}" : ''
  """
  wget -O records.embl '${params.databases.ena_coding.url}/${source.accessions}?download=true'
  rnac ena parse-coding --taxid ${source.taxid} $assembly $end records.embl .
  """
}

workflow ena_coding {
  main:
    channel.fromList(params.databases.ena_coding.sources) | parse_ena_coding | set { data }
  emit: data
}
