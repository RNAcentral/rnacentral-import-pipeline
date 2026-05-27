process fetch_model_stats {
  container params.r2dt.container
  containerOptions "${params.common_container}"

  input:
  val(_flag)

  output:
  path('info.csv'), emit: info
  path '*.tsv', emit: metadata

  script:
  """
  find /rna/r2dt/data -type f -name '*.cm' | xargs -I {} cmstat {}  | grep -v ^\\# | awk '{ printf("%s,%d,%d\n", \$3, \$6, \$8); }' | sort -u > info.csv
  cp /rna/r2dt/data/rnasep/metadata.tsv rnase-p.tsv
  cp /rna/r2dt/data/crw-metadata.tsv crw.tsv
  awk 'FNR==1 && NR>1{next}1' /rna/r2dt/data/ribovision*/metadata.tsv > ribovision.tsv
  """
}

process create_model_info {
  input:
  tuple val(model_source), path(metadata), path(info)

  output:
  path('models.csv')

  script:
  def args = model_source == "crw" ? "$info $metadata" : "$metadata"
  """
  rnac r2dt model-info $model_source $args models.csv
  """
}

process store_model_info {
  input:
  path('models*.csv')
  path(load)

  output:
    val('model info stored')

  script:
  """
  pgloader --on-error-stop $load
  """
}

workflow model_info {
  take: ready
  main:
    Channel.fromPath('files/r2dt/load-models.ctl') | set { load }

    fetch_model_stats(ready)

    fetch_model_stats.out.metadata \
    | flatten \
    | map { fn -> [fn.baseName, fn] } \
    | combine(fetch_model_stats.out.info) \
    | create_model_info \
    | collect \
    | set { model_info }

    store_model_info(model_info, load) | set { done }
  emit: done
}

workflow {
  model_info(Channel.of('ready'))
}
