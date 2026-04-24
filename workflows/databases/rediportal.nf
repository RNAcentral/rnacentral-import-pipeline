process build_rnc_bedfile {
  tag { assembly_id }

  input:
    tuple val(assembly_id), path(query)

  output:
    tuple val(assembly_id), path("rnc_regions.bed")

  """
  set -euo pipefail
  psql -v ON_ERROR_STOP=1 -v "assembly_id=$assembly_id" -f $query "$PGDATABASE" > result.json

  rnac ftp-export coordinates as-bed result.json |\
  sort -k1,1 -k2,2n > rnc_regions_gene.bed
  bed-expander rnc_regions_gene.bed rnc_regions.bed
  """
}

process fetch_rediportal_inputs {
  tag { assembly_id }

  input:
    tuple val(assembly_id), val(genome_build), val(remote)

  output:
    tuple val(assembly_id), val(genome_build), path("rediportal.txt"), path("rediportal.bed")

  """
  set -euo pipefail
  curl -L "$remote" | gzip -d > rediportal.txt
  rnac rediportal build-bed --genome-build "$genome_build" rediportal.txt rediportal.bed
  """
}

process intersect_rnc_rediportal {
  tag { assembly_id }

  memory '16GB'

  input:
    tuple val(assembly_id), path(rnc_bed), val(genome_build), path(redi_meta), path(redi_bed)

  output:
    path("features.csv")

  """
  rnac rediportal parse-bed --genome-build "$genome_build" $redi_bed $redi_meta $rnc_bed features.csv
  """
}

process load_rediportal {

  memory 6.GB

  input:
    tuple path(features), path(ctl)

  output:
    val('rediportal done')

  """
  split-and-load $ctl *.csv ${params.databases.rediportal.chunk_size} rediportal-data
  """

}

workflow rediportal {
  take: ready
  emit: done
  main:
    if( params.databases.rediportal.run ) {
      Channel.fromPath("files/ftp-export/genome_coordinates/query.sql") | set { region_query }
      Channel.fromPath("files/rediportal/load.ctl") | set { load_query }
      Channel
        .fromList(params.databases.rediportal.inputs)
        .map { input -> tuple(input.assembly_id, input.genome_build, input.remote) }
        | set { redi_inputs }

      redi_inputs
        .map { assembly_id, _genome_build, _remote -> assembly_id }
        .unique()
        .combine(region_query)
        | build_rnc_bedfile
        | set { rnc_bedfile }

      redi_inputs | fetch_rediportal_inputs | set { redi_data }

      rnc_bedfile
        .combine(redi_data)
        .filter { left_assembly, _rnc_bed, right_assembly, _genome_build, _redi_meta, _redi_bed ->
          left_assembly == right_assembly
        }
        .map { _left_assembly, rnc_bed, assembly_id, genome_build, redi_meta, redi_bed ->
          tuple(assembly_id, rnc_bed, genome_build, redi_meta, redi_bed)
        }
        | intersect_rnc_rediportal
        | combine(load_query)
        | load_rediportal
        | set { done }

    }
    else {
      Channel.of('rediportal not run') | set { done }
    }

}


workflow {
  rediportal(Channel.of('ready'))
}
