#!/usr/bin/env nextflow

// process find_chunks {
//   output:
//   stdout raw_ranges

//   """
//   rnac upi-ranges ${params.precompute.max_entries}
//   """
// }

// raw_ranges
//   .splitCsv()
//   .combine(Channel.fromPath('files/precompute/query.sql'))
//   .set { ranges }

// process precompute_range_query {
//   beforeScript 'slack db-work precompute-range || true'
//   afterScript 'slack db-done precompute-range || true'
//   maxForks params.precompute.maxForks

//   input:
//   set val(min), val(max), file(query) from ranges

//   output:
//   file 'raw-precompute.json' into precompute_raw

//   """
//   psql --variable min=$min --variable max=$max -f "$query" '$PGDATABASE' > raw-precompute.json
//   """
// }

// process precompute_range {
//   memory '4 GB'

//   input:
//   file(raw) from precompute_raw

//   output:
//   file 'precompute.csv' into precompute_results
//   file 'qa.csv' into qa_results

//   """
//   rnac precompute from-file $raw
//   """
// }

// process load_precomputed_data {
//   beforeScript 'slack db-work loading-precompute || true'
//   afterScript 'slack db-done loading-precompute || true'

//   input:
//   file('precompute*.csv') from precompute_results.collect()
//   file('qa*.csv') from qa_results.collect()
//   file pre_ctl from Channel.fromPath('files/precompute/load.ctl')
//   file qa_ctl from Channel.fromPath('files/precompute/qa.ctl')
//   file post from Channel.fromPath('files/precompute/post-load.sql')

//   // Use cp to copy the ctl files and not stageInMode because we don't want to
//   // have to copy a large number of large files before running this. We only
//   // need a copy of these two files in the local directory. This is because
//   // pgloader will look at the location of the end of symlink when trying to
//   // find files and not the current directory or where the symlink lives.
//   """
//   cp $pre_ctl _pre.ctl
//   cp $qa_ctl _qa.ctl
//   pgloader _pre.ctl
//   pgloader _qa.ctl
//   psql -f $post "$PGDATABASE"
//   """
// }

process mods_for_feedback {
  input:
  file(query) from Channel.fromPath('files/precompute/find-mod-info.sql')

  output:
  file(info) into raw_mods

  """
  psql -f "$query" "$PGDATABASE" > info
  """
}

raw_mods
  .splitCsv()
  .combine(Channel.fromPath('files/ftp-export/genome_coordinates/query.sql'))
  .set { mods }

process generate_feedback_report {
  input:
  set val(assembly), val(mod), file(query) from mods

  output:
  file('combined.tsv') into feedback

  shell:
  genome="${mod}.bed"
  '''
  psql -v "assembly_id=!{assembly}" -f !{query} "$PGDATABASE" > result.json
  rnac ftp-export coordinates as-bed result.json !{genome}
  cut -f15 !{genome} | tr ',' '\\n' | sort -u | grep -vi !{mod} > db_list.txt

  for db in `cat db_list.txt`; do
      grep $db !{genome} > ${db}.bed
      bedtools intersect -wa -wb -a ${db}.bed -b !{genome} > overlap-${db}.bed
      bedtools subtract -A -a !{genome} -b ${db}.bed > no-overlap-${db}.bed
      db_report=${db}-report.tsv

      awk '{split($30, dbs, ","); for(i in dbs) print $4, "overlap", dbs[i];}' overlap-${db}.bed | sort -u >> $db_report
      awk -v awk_db=$db '{print $4, "no_overlap", awk_db}' no-overlap-${db}.bed | sort -u >> $db_report
  done

  sort -u *-report.tsv > combined.tsv
  '''
}

process import_feedback {
  echo true

  input:
  file('feedback*.tsv') from feedback.collect()
  file(ctl) from Channel.fromPath('files/precompute/feedback.ctl')

  """
  cp $ctl local_$ctl
  cat feedback*.tsv > merged.tsv
  pgloader --on-error-stop local_$ctl
  """
}
