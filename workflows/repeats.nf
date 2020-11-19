process find_genomes_with_repeats {
  when { params.precompute.run }

  input:
  tuple path(query), path(connections)

  output:
  path("info.csv")

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > assemblies.csv
  rnac repeats find-databases $connections assemblies.csv info.csv
  """
}

process query_ensembl {
  maxForks 10

  input:
  tuple val(assembly), val(conn_name), val(database), path(query)

  output:
  tuple val(assembly), path('raw.bed')

  script:
  def conn = params.connections[conn_name]
  """
  mysql -N --host $conn.host --port $conn.port --user $conn.user --database $database < $query > raw.bed
  """

}

process fetch_ensembl_data {
  tag { "$species-$assembly" }
  memory '4GB'

  input:
  tuple val(assembly), path(raw)

  output:
  path("$assembly-repeats"), emit: repeats

  script:
  def out = "$assembly-repeats"
  def coord = "${assembly}.bed.bgz"
  """
  mkdir $out
  pushd $out
  bedtools merge -i ../$raw | sort -k1V -k2n -k3n | bgzip > $coord
  tabix -s 1 -b 2 -e 3 $coord
  rnac repeats build-info-directory --chromosome-column 1 --start-column -2 --stop-column 3 $assembly .
  popd
  """
}

workflow repeats {
  emit: repeats
  main:
    Channel.fromPath('files/repeats/find-assemblies.sql') \
    | combine(Channel.fromPath('config/databases.json')) \
    | find_genomes_with_repeats \
    | splitCsv \
    | combine(Channel.fromPath('files/repeats/extract-repeats.sql')) \
    | query_ensembl \
    | fetch_ensembl_data

    fetch_ensembl_data.out.repeats | collect | set { repeats }
}
