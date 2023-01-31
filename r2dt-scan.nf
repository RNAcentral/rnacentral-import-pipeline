#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { r2dt } from './workflows/r2dt'

process extract_gtrnadb_metadata {
  container params.r2dt.container
  memory '256 MB'
  input:
    path(model_path)
  output:
    path("basepairs.csv")


  shell:
  '''
  cmstat !{model_path} | awk ' /^[^#]/ { sep=","; printf "%s%s%s\\n", $2,sep,$8 }' > basepairs.csv
  '''
}

process parse_gtrnadb_model {
  memory '256 MB'

  input:
    tuple path(model_path), path(basepairs)
  output:
    path("model_data.csv")

  shell:
    '''
    sort -k 1 !{basepairs} > basepairs_sorted.csv
    rnac r2dt model-info gtrnadb !{model_path} model_data_nbp.csv
    join -t","  model_data_nbp.csv basepairs_sorted.csv -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2 > model_data.csv
    '''
}

process extract_ribovision_metadata {
    container params.r2dt.container
    memory '256 MB'

    output:
      path('length_basepair.csv')

    shell:
    '''
      cmstat /rna/r2dt/data/ribovision-lsu/cms/all.cm | awk '/^[^#]/ {sep=","; printf "%s%s%s%s%s\\n",$2,sep,$6,sep,$8}' > basepair_length_lsu
      cmstat /rna/r2dt/data/ribovision-ssu/cms/all.cm | awk '/^[^#]/ {sep=","; printf "%s%s%s%s%s\\n",$2,sep,$6,sep,$8}' > basepair_length_ssu
      cat basepair_length_lsu basepair_length_ssu  | sort -k 1 > length_basepair.csv
  '''
}

process parse_ribovision_models {
  memory '256 MB'

  input:
    tuple val(ribovision_metadata_url), path(length_basepair)



  output:
    path("model_data.csv")

  shell:
  '''
  wget !{ribovision_metadata_url}
  rnac r2dt model-info ribovision metadata.tsv model_data_us.csv
  sort -k 1 model_data_us.csv > model_data_s.csv
  join -t","  model_data_s.csv !{length_basepair} -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3 > model_data.csv
  '''

}

process extract_rnasep_metadata {
  container params.r2dt.container
  memory '256 MB'

  output:
    path('length_basepair.csv')

  shell:
  '''
  cmstat /rna/r2dt/data/rnasep/cms/all.cm | awk '/^[^#]/ {sep=","; printf "%s%s%s%s%s\\n",$2,sep,$6,sep,$8}' > length_basepair.csv
  '''
}

process parse_rnasep_models {
  memory '256 MB'

  input:
    tuple val(rnasep_metadata_url), path(length_basepair)

  output:
    path("model_data.csv")

  shell:
  '''
  wget !{rnasep_metadata_url}
  sed -i 's/\\tNRC-1\\t/\\t/g' metadata.tsv
  rnac r2dt model-info rnase-p metadata.tsv model_data_us.csv
  sort -k1 model_data_us.csv > model_data_s.csv
  join -t","  model_data_s.csv !{length_basepair} -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3 > model_data.csv
  '''

}

process extract_rfam_metadata {
  container params.r2dt.container
  memory '256 MB'

  input:
    path(all_models)

  output:
    path('basepairs.csv')

  shell:
  '''
  cmstat !{all_models} | awk '/^[^#]/ {sep=","; printf "%s%s%s\\n",$3,sep,$8}' > basepairs.csv
  '''
}

process parse_rfam_models {
  memory '256 MB'

  input:
    tuple path(all_models), path(basepairs)
  output:
    path("model_data.csv")

  shell:
  '''
  rnac r2dt model-info rfam !{all_models} $PGDATABASE model_data_nbp.csv
  join -t","  model_data_nbp.csv !{basepairs} -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2 > model_data.csv
  '''
}

process extract_crw_metadata {
  container params.r2dt.container
  memory '256 MB'

  input:
    path(all_models)

  output:
    path("basepairs.csv")

  shell:
  '''
    cmstat !{all_models} | awk '/^[^#]/ {sep=","; printf "%s%s%s\\n",$2,sep,$8}' > basepairs.csv
  '''
}

process parse_crw_models {
  memory '256 MB'

  input:
    tuple path(all_models), val(metadata), path(basepairs)
  output:
    path("model_data.csv")

  script:
  """
  wget $metadata -O metadata.tsv
  sed -i 's/taxid  rna_type/taxid\trna_type/g' metadata.tsv
  rnac r2dt model-info crw $all_models metadata.tsv model_data_nbp.csv
  join -t","  model_data_nbp.csv $basepairs -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2 > model_data.csv
  """
}

process load_models {

  memory '256 MB'

  input:
    path(all_data)
    path(ctl)

  output:
    val('models loaded')

  script:
  """
  split-and-load $ctl $all_data ${params.import_data.chunk_size}
  """
}





workflow {
  rfam_models = Channel.of("$baseDir/singularity/bind/r2dt/data/cms/rfam/all.cm")
  crw_models = Channel.of("$baseDir/singularity/bind/r2dt/data/cms/crw/all.cm")
  crw_metadata = Channel.of("https://raw.githubusercontent.com/RNAcentral/R2DT/v1.3/data/crw-metadata.tsv")
  gtrnadb_models = Channel.fromPath("$baseDir/singularity/bind/r2dt/data/cms/gtrnadb/*.cm")
  ribovision_lsu_metadata_url = Channel.of("https://raw.githubusercontent.com/RNAcentral/R2DT/v1.3/data/ribovision-lsu/metadata.tsv")
  ribovision_ssu_metadata_url = Channel.of("https://raw.githubusercontent.com/RNAcentral/R2DT/v1.3/data/ribovision-ssu/metadata.tsv")

  rnasep_metadata_url = Channel.of("https://raw.githubusercontent.com/RNAcentral/R2DT/v1.3/data/rnasep/metadata.tsv")

  load_ctl = Channel.of("$baseDir/files/r2dt/load-models.ctl")

  rfam_models | extract_rfam_metadata | set { rfam_basepairs }
  rfam_models.combine(rfam_basepairs) | parse_rfam_models | set { rfam_data }

  crw_models | extract_crw_metadata | set { crw_basepairs }
  crw_models.combine(crw_metadata) | combine(crw_basepairs) | parse_crw_models | set { crw_data }

  gtrnadb_models | extract_gtrnadb_metadata | collectFile() {csvfile -> [csvfile.name, csvfile.text]} | set { gtrnadb_basepairs }
  gtrnadb_models.combine(gtrnadb_basepairs) | parse_gtrnadb_model | collectFile() {csvfile -> [csvfile.name, csvfile.text]} | set { gtrnadb_data }

  extract_ribovision_metadata | set { ribovision_length_basepair }
  ribovision_lsu_metadata_url.mix(ribovision_ssu_metadata_url) | combine(ribovision_length_basepair) | parse_ribovision_models | set { ribovision_data }

  extract_rnasep_metadata | set { rnasep_length_basepair }
  rnasep_metadata_url.combine(rnasep_length_basepair) | parse_rnasep_models | set { rnasep_data }


  rfam_data.mix(crw_data, gtrnadb_data, ribovision_data, rnasep_data) | collectFile() {csvfile -> [csvfile.name, csvfile.text]} | set { all_data }


  load_models(all_data, load_ctl) | set { model_load }

  model_load | r2dt | set { done }

}
