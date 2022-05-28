process generate_files {
  when { params.rfam.run }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(_flag)

  output:
  path 'rfam', emit: cm_files
  path 'version_file', emit: version_info

  script:
  def base = params.rfam.files
  """
  mkdir rfam
  cd rfam
  wget -O Rfam.clanin "$base/Rfam.clanin"
  wget -O Rfam.cm.gz "$base/Rfam.cm.gz"
  gzip -d *.gz
  cmpress Rfam.cm
  cd ..
  wget -O version_file "$base/README*"
  """
}

process sequences {
  memory '20GB'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  // clusterOptions '-R "rusage[scratch=4000]"'

  input:
  tuple path(version), path(active_xrefs), path(computed), path(compute_missing)

  output:
  tuple path(version), path('parts/*.fasta')

  """
  export TMPDIR="$baseDir/work/tmp"
  psql -v ON_ERROR_STOP=1 -f "$active_xrefs" "$PGDATABASE" | sort -u > active-urs
  psql -v ON_ERROR_STOP=1 -v 'name=rfam' -f "$computed" "$PGDATABASE" | sort > computed
  comm -23 active-urs computed > urs-to-compute
  psql -q -v ON_ERROR_STOP=1 -f "$compute_missing" "$PGDATABASE" > raw.json
  mkdir parts
  split --filter 'json2fasta --only-valid-easel - - >> \$FILE.fasta' --lines ${params.rfam.chunk_size} raw.json parts/
  """
}

process scan {
  cpus { params.rfam.cpus }
  memory { params.rfam.memory * params.rfam.cpus }
  errorStrategy 'ignore'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple path(version), path('sequences.fasta'), path(cm_files)

  output:
  path 'hits.csv', emit: hits
  path 'attempted.csv', emit: attempted

  """
  cmscan \
    -o output.inf \
    --tblout results.tblout \
    --clanin $cm_files/Rfam.clanin \
    --oclan \
    --fmt 2 \
    --acc \
    --cut_ga \
    --rfam \
    --notextw \
    --nohmmonly \
    "$cm_files/Rfam.cm" \
    sequences.fasta

  rnac qa rfam results.tblout hits.csv
  rnac qa create-attempted sequences.fasta rfam $version attempted.csv
  """
}

process import_data {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('raw*.csv')
  path(ctl)
  path('attempted*.csv')
  path('attempted.ctl')

  output:
  val('rfam done')

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} rfam
  split-and-load attempted.ctl 'attempted*.csv' ${params.import_data.chunk_size} attempted-rfam
  """
}

workflow rfam_scan {
  take: ready
  emit: done
  main:
    Channel.fromPath("files/find-active-xrefs-urs.sql").set { active_xref_sql }
    Channel.fromPath("files/qa/computed.sql").set { computed_sql }
    Channel.fromPath("files/qa/compute-required.sql").set { compute_required_sql }
    Channel.fromPath("files/rfam-scan/load.ctl").set { ctl }
    Channel.fromPath("files/rfam-scan/load-attempted.ctl").set { attempted_ctl }

    generate_files(ready)

    generate_files.out.version_info \
    | combine(active_xref_sql) \
    | combine(computed_sql) \
    | combine(compute_required_sql) \
    | sequences \
    | flatMap { version, files ->
      (files instanceof ArrayList) ? files.collect { [version, it] } : [[version, files]]
    } \
    | filter { v, f -> !f.isEmpty() } \
    | combine(generate_files.out.cm_files) \
    | scan

    scan.out.hits | collect | set { hits }
    scan.out.attempted | collect | set { attempted }

    import_data(hits, ctl, attempted, attempted_ctl) | set { done }
}
