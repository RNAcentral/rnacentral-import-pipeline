process fetch_and_process {
  output:
  path('*.{csv,parquet}')

  when: params.databases.ribovision?.run

  script:
  """
  git clone --depth 1 --filter=blob:none --sparse "$params.databases.ribovision.r2dt_repo" r2dt
  git -C r2dt sparse-checkout set data/ribovision-lsu data/ribovision-ssu

  rnac ribovision r2dt-to-fasta r2dt/data sequences.fasta
  rnac ribovision parse r2dt/data sequences.fasta
  """
}

workflow ribovision {
  main:
    fetch_and_process | set { data }
  emit: data
}
