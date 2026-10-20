process fetch_species {
  tag { "$key" }

  input:
  tuple val(key), val(bed)

  output:
  tuple val(key), path("${key}.fa"), path("${key}.bed")

  when: params.databases.noncode?.run

  script:
  def base = params.databases.noncode.remote
  // Every species has sequences; only human and mouse have coordinates, so the
  // others get an empty BED to keep the tuple shape.
  def fetch_bed = bed ? "wget -O ${key}.bed.gz '$base/${bed}' && gzip -df ${key}.bed.gz" : "touch ${key}.bed"
  """
  wget -O ${key}.fa.gz '$base/NONCODEv6_${key}.fa.gz'
  gzip -df ${key}.fa.gz
  $fetch_bed
  """
}

process parse_species {
  tag { "$key" }
  memory '4GB'

  input:
  tuple val(key), path(fasta), path(bed)

  output:
  path('*.{csv,parquet}')

  script:
  def with_bed = bed.size() > 0 ? "--bed $bed" : ""
  """
  rnac noncode parse $with_bed $key $fasta .
  """
}

workflow noncode {
  main:
    channel.fromList(params.databases.noncode.species) \
    | fetch_species \
    | parse_species \
    | set { data }

  emit: data
}
