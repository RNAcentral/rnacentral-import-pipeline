process fetch_model_stats {
  when { params.r2dt.run }
  container params.r2dt.container
  containerOptions "--bind ${params.r2dt.cms_path}:/rna/r2dt/data/cms" 

  input:
  val(_flag)

  output:
  path('info.csv'), emit: info
  path '*.tsv', emit: metadata

  """
  find /rna/r2dt/data/cms -type f -name '*.cm' | xargs -I {} cmstat {}  | grep -v ^\\# | awk '{ printf("%s,%d,%d\n", \$3, \$6, \$8); }' | sort -u > info.csv
  cp /rna/r2dt/data/rnasep/metadata.tsv rnasep.tsv
  cp /rna/r2dt/data/crw-metadata.tsv crw.tsv
  cat /rna/r2dt/data/ribovision*/metadata.tsv > ribovision.tsv
  """
}

process create_model_info {
  input:
  tuple val(model_source), path(metadata), path(info)

  """
  rnc r2dt model-info $model_source $info $metadata models.csv
  """
}

process store_model_info {
  input:
  path('models*.csv')
  path(load)

  """
  pgloader --on-error-stop $load
  """
}

workflow model_info {
  take: ready
  emit: done
  main:
    Channel.fromPath('files/r2dt/load-models.ctl') | set { load }

    fetch_model_stats(ready)

    fetch_model_stats.out.metadata \
    | map { fn -> [fn.baseName, fn] } \
    | combine(fetch_model_stats.out.info) \
    | create_model_info \
    | collect \
    | set { model_info }

    store_model_info(model_info, load)
}
