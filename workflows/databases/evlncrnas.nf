process fetch {
  when { params.databases.evlncrnas.run }

  output:
  path('EVLncRNAs2_alldata')

  """
  wget --no-check-certificate --read-timeout=30 -t 1 http://www.sdklab-biophysics-dzu.net/EVLncRNAs2/Public/home/download/EVLncRNAs2_alldata.zip
  unzip EVLncRNAs2_alldata.zip
  """

}

process rnc_dump {
  when { params.databases.evlncrnas.run }

  input:
  path(query)

  output:
  path('*.csv')

  """
  psql -f $query $PGDATABASE > ev_lookup.csv
  """
}


process parse {
  memory '16 GB'

  input:
  tuple path(ev_data), path(rnc_data)

  output:
  path('*.csv')

  """
  rnac evlncrnas parse $ev_data $rnc_data .
  """
}


workflow evlncrnas {
  emit: data
  main:
  Channel.fromPath('files/import-data/evlncrnas/dump-lookup.sql') | set { dump_sql }
    fetch | set {ev_data}
    dump_sql | rnc_dump | set {rnc_data}

    ev_data.combine(rnc_data) | parse | set {data}
}
