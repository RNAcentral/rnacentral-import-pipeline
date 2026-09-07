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
  // The 144 SSU file is 3.5Gb gzipped and OOMed on the 1Gb cluster default.
  // Scale on attempt rather than exit status: a SLURM kill often reports none.
  memory { 8.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 2

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
