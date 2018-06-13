#!/usr/bin/env nextflow

sequences_query = Channel.from('files/qa/rfam-scan.sql')
process fetch_sequences {
  input:
  file query from sequences_query

  output:
  file 'rnacentral.fasta' into sequences_to_split

  """
  psql -f $query | tsv2fasta.py - rnacentral.fasta
  """
}

process randomize_sequences {
  input:
  file sequences from sequences_to_split

  output:
  file('**/rnac.fa') into sequences_to_scan

  """
  esl-randomize-sqfile.pl -O rnac.fa -L -N 3000 $sequences 1.0
  """
}

cm_files = Channel.fromPath(params.qa.rfam_scan.cm_file)
process infernal_scan {
  queue 'mpi-rh7'
  cpus params.qa.rfam_scan.cpus
  clusterOptions "-M ${params.qa.rfam_scan.cm_memory} -R 'rusage[mem=${params.qa.rfam_scan.cm_memory}]' -R 'span[hosts=1]' -a mpi mpiexec -mca btl ^openbib -np ${params.rfam_anotations.cpus}"

  input:
  file sequences from sequences_to_scan
  file cm_file from cm_files

  output:
  file 'results.tblout' into infernal_results

  """
  cmscan \
    -o output.inf \
    --tblout results.tblout \
    --clanin ${params.qa.rfam_scan.clans}
    --oclan \
    --fmt 2 \
    --acc \
    --cut_ga \
    --rfam \
    --notextw \
    --nohmmonly \
    --mpi \
    $cm_file \
    $sequences
  """
}

process process_hits {
  input:
  file('results.tblout*') from infernal_results.collect()

  output:
  file('hits.csv') into processed_hits

  """
  cat results.tblout* > all_results.tblout
  rnac rfam-anntotations process all_result.tblout hits.csv
  """
}

hits_ctl = Channel.from('files/qa/rfam-scan.ctl')
process import_hits {
  input:
  file('hits.csv') from processed_hits
  file hit_ctl from hits_cl

  """
  pgloader $hit_ctl
  """
}
