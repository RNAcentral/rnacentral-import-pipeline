process fetch {
  output:
  path('*.rnac')

  when: params.databases.silva?.run

  script:
  """
  # before SILVA 144 split SSU and LSU across releases:
  # wget -e robots=off -nH -r --cut-dirs 3 --no-parent -A "SILVA_*Parc.rnac.gz" \$params.databases.silva.remote
  wget ${params.databases.silva.remote.join(' ')}
  gzip -d *.gz
  """
}

process parse {
  tag { "$raw.name" }

  input:
  tuple path(raw), path(taxonomy)

  output:
  path('*.{csv,parquet}')

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
