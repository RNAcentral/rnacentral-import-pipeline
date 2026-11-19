process parse_refseq_coding {
  tag { source.name }
  memory '4GB'

  input:
  val source

  output:
  path('*.{csv,parquet}')

  when: params.databases.refseq_coding?.run

  script:
  def config = params.databases.refseq_coding
  def url = source.accessions ? "${config.efetch}&id=${source.accessions}" : "${config.genomes}/${source.path}"
  def assembly = source.assembly ? "--assembly '${source.assembly}'" : ''
  """
  wget -O records '$url'
  gzip -dcf records > records.gbff
  rnac refseq parse-coding --taxid ${source.taxid} $assembly records.gbff .
  """
}

workflow refseq_coding {
  main:
    channel.fromList(params.databases.refseq_coding.sources) | parse_refseq_coding | set { data }
  emit: data
}
