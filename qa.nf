process fetch_rfam_files {
  input:
  val(base)

  output:
  path(name), emit: scan_files
  path('version_file'), emit: rfam_version

  script:
  """
  mkdir rfam
  pushd rfam
  fetch generic "$base/Rfam.clanin" Rfam.clanin
  fetch generic "$base/Rfam.cm.gz" Rfam.cm.gz
  gzip -d *.gz
  cmpress Rfam.cm
  popd

  fetch generic "$base/README" version_file
  """
}

process fetch_rfam_sequences {
  memory '5GB'

  input:
  tuple path(version), path(active_xrefs), path(computed), path(compute_missing)

  output:
  tuple path(version), path('parts/*.fasta') into split_qa_sequences

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$active_xrefs" "$PGDATABASE" | sort -u > active-urs
  psql -v ON_ERROR_STOP=1 -v 'name=rfam' -f "$computed" "$PGDATABASE" | sort > computed
  comm -23 active-urs computed > urs-to-compute
  psql -q -v ON_ERROR_STOP=1 -f "$compute_missing" "$PGDATABASE" > raw.json
  mkdir parts
  split --filter 'json2fasta.py --only-valid-easel - \$FILE.fasta' --lines ${params.qa[name].chunk_size} raw.json parts/
  """
}

process rfam_scan {
  cpus { params.qa.rfam.cpus }
  memory { params.qa.rfam.memory * params.qa.rfam.cpus }
  errorStrategy 'ignore'

  input:
  tuple path(version), path('sequences.fasta'), path(dir)

  output:
  path('hits.csv'), emit: hits
  path('attempted.csv'), emit: attempted

  script:
  """
  cmscan \
    -o output.inf \
    --tblout results.tblout \
    --clanin $dir/Rfam.clanin \
    --oclan \
    --fmt 2 \
    --acc \
    --cut_ga \
    --rfam \
    --notextw \
    --nohmmonly \
    "$dir/Rfam.cm" \
    sequences.fasta
  rnac qa rfam results.tblout hits.csv

  rnac qa create-attempted sequences.fasta rfam $version attempted.csv
  """
}

process import_rfam_data {
  input:
  path('raw*.csv')
  path(ctl)
  path('attempted*.csv')
  path('attempted.ctl')

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  split-and-load attempted.ctl 'attempted*.csv' ${params.import_data.chunk_size} attempted-$name
  """
}

workflow rfam {
  Channel.fromPath('files/qa/rfam.ctl').set { hits_ctl }
  Channel.fromPath('files/qa/attempted/rfam.ctl').set { attempted_ctl }

  fetch_rfam_files(params.qa.rfam.files)

  fetch_rfam_files.out.rfam_version \
  | combine(Channel.fromPath("files/find-active-xref-urs.sql")) \
  | combine(Channel.fromPath("files/qa/computed.sql")) \
  | combine(Channel.fromPath("files/qa/compute-required.sql")) \
  | fetch_rfam_sequences \
  | join(fetch_rfam_files.out.scan_files) \
  | flatMap { version, fns, extra -> fns.inject([]) { acc, fn -> acc << [version, fn, extra] } } \
  | rfam_scan

  rfam_scan.out.hits | collect | set { hits }
  rfam_scan.out.attempted | collect | set { attempted}

  import_rfam_data(hits, hits_ctl, attempted, attempted_ctl)
}

workflow qa {
  rfam()
}

workflow {
  qa()
}
