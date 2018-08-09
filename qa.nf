#!/usr/bin/env nextflow

process fetch_sequences {
  input:
  file query from Channel.fromPath('files/qa/rfam-scan.sql')

  output:
  file 'rnacentral.fasta' into sequences_to_split

  """
  psql -f "$query" "$PGDATABASE" | json2fasta.py - rnacentral.fasta
  """
}

process randomize_sequences {
  input:
  file sequences from sequences_to_split

  output:
  file('rnac_*.fa') into chunks

  """
  esl-randomize-sqfile.pl -O rnac.fa -L -N ${params.qa.rfam_scan.chunk_size} $sequences 1.0
  """
}

sequences_to_scan = chunks.flatMap()
cm_files = Channel.fromPath(params.qa.rfam_scan.cm_files)
process infernal_scan {
  queue 'mpi-rh7'
  cpus params.qa.rfam_scan.cpus
  clusterOptions "-M ${params.qa.rfam_scan.cm_memory} -R 'rusage[mem=${params.qa.rfam_scan.cm_memory}]' -a openmpi"
  module 'mpi/openmpi-x86_64'

  input:
  file 'sequences.fasta' from sequences_to_scan
  file cm_file from cm_files.collect()

  output:
  file 'results.tblout' into infernal_results

  """
  mpiexec -mca btl ^openbib -np ${params.qa.rfam_scan.cpus} \
  cmscan \
    -o output.inf \
    --tblout results.tblout \
    --clanin ${params.qa.rfam_scan.clans} \
    --oclan \
    --fmt 2 \
    --acc \
    --cut_ga \
    --rfam \
    --notextw \
    --nohmmonly \
    --mpi \
    "Rfam.cm" \
    sequences.fasta
  """
}

process process_hits {
  input:
  file('results.tblout*') from infernal_results.collect()

  output:
  file('hits.csv') into processed_hits

  """
  cat results.tblout* | rnac qa tblout2csv - hits.csv
  """
}

hits_ctl = Channel.fromPath('files/qa/rfam-scan.ctl')
process import_hits {
  input:
  file('hits.csv') from processed_hits
  file hit_ctl from hits_ctl

  """
  pgloader $hit_ctl
  """
}
