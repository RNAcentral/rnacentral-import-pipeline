process fetch {
  output:
  path('EVLncRNAs3_alldata')

  when: params.databases.evlncrnas?.run

  script:
  """
  wget --no-check-certificate --read-timeout=30 -t 1 \
    https://request.sdklab-biophysics-dzu.net/uploads/files/EVLncRNAs3_alldata.zip
  mkdir -p EVLncRNAs3_alldata
  unzip EVLncRNAs3_alldata.zip -d EVLncRNAs3_alldata
  """

}

process rnc_dump {
  input:
  path(query)

  output:
  path('*.{csv,parquet}')

  when: params.databases.evlncrnas?.run

  script:
  """
  psql -f $query \$PGDATABASE > ev_lookup.csv
  """
}


process map_diseases {
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple path(ev_data), path(overrides)

  output:
  path('disease_mapping.tsv')

  script:
  """
  rnac evlncrnas disease-terms $ev_data disease_terms.txt
  rnac ols map-text --overrides $overrides disease_terms.txt disease_mapping.tsv
  """
}


process parse {
  memory '16 GB'

  input:
  tuple path(ev_data), path(rnc_data), path(disease_mapping)

  output:
  path('*.{csv,parquet}')

  script:
  """
  rnac evlncrnas parse --disease-mapping $disease_mapping $ev_data $rnc_data .
  """
}


workflow evlncrnas {
  main:
  channel.fromPath('files/import-data/evlncrnas/dump-lookup.sql') | set { dump_sql }
    fetch | set {ev_data}
    dump_sql | rnc_dump | set {rnc_data}
    channel.fromPath('files/disease-mapping-overrides.tsv') | set { overrides }
    ev_data.combine(overrides) | map_diseases | set {disease_mapping}

    ev_data.combine(rnc_data).combine(disease_mapping) | parse | set {data}
  emit: data
}
