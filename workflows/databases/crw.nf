process fetch_and_process {
  output:
  path('*.{csv,parquet}')

  when: params.databases.crw?.run

  script:
  """
  git clone --depth 1 --filter=blob:none --sparse "$params.databases.crw.r2dt_repo" r2dt
  git -C r2dt sparse-checkout set data/crw-fasta data/crw-metadata.tsv

  rnac crw r2dt-to-fasta r2dt/data/crw-fasta sequences.fasta
  rnac crw parse r2dt/data/crw-metadata.tsv sequences.fasta
  """
}

workflow crw {
  main:
    fetch_and_process | set { data }
  emit: data
}
