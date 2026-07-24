process fetch {
  output:
  path('*.rnac')

  when: params.databases.silva?.run

  script:
  """
  wget -e robots=off -nH -r --cut-dirs 3 --no-parent -A "SILVA_*Parc.rnac.gz" $params.databases.silva.remote
  gzip -d *.gz
  """
}

process parse {
  tag { "$raw.name" }

  input:
  tuple path(raw), path(taxonomy)

  output:
  path('*.csv')

  script:
  """
  rnac silva parse $raw $taxonomy .
  """
}

workflow silva {
  take: tax_info
  main:
    fetch \
    | flatten \
    | combine(tax_info) \
    | parse \
    | set { data }
  emit: data
}
