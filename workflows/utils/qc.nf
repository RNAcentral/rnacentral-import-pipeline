nextflow.enable.dsl=2

// Databases with run=true this run (mirrors main.nf; ensembl special-cased),
// normalised for matching rnc_database.descr/display_name in SQL.
def running_databases() {
  def out = []
  def dbs = params.databases ?: [:]
  def ensembl_descr = [
    vertebrates: 'ensembl', plants: 'ensembl_plants', fungi: 'ensembl_fungi',
    protists: 'ensembl_protists', metazoa: 'ensembl_metazoa',
  ]
  dbs.each { key, db ->
    if (!(db instanceof Map)) return
    if (key == 'ensembl') {
      ensembl_descr.each { sub, descr -> if (db[sub]?.get('run', false)) out << descr }
    } else if (db.get('run', false)) {
      out << key
    }
  }
  return out
}

def run_dbs_arg() {
  running_databases().collect { it.toLowerCase().replaceAll('[_ ]', '') }.unique().join(',')
}

def qc_dir() { params.qc?.publish_dir ?: 'qc-reports' }

// Echo a finished report to the Nextflow console (publishDir handles the file).
def print_report(report) {
  def dir = qc_dir()
  report.view { f -> "\n— QC report — also written to ${dir}/${f.name}\n${f.text}" }
}

// Shared report header for every stage: title, UTC timestamp, target database
// (profile label pro/dev/tst + the live host:port/db it actually connected to,
// so a wrong-DB mix-up is visible), and wall-clock elapsed since pipeline start.
// Returns bash that sets \$target/\$elapsed then echoes the lines — drop the
// interpolation inside the report's { ... } block.
def qc_header(stage) {
  def env_label = (workflow.profile ?: '').tokenize(',')*.trim().find { it in ['pro', 'dev', 'tst'] } ?: ''
  def prefix = env_label ? env_label + ' — ' : ''
  def run_start = workflow.start ?: ''
  return """\
    target=\$(psql -Aqt -c "select coalesce(host(inet_server_addr()),'localhost')||':'||coalesce(inet_server_port()::text,'?')||'/'||current_database()" "\$PGDATABASE" 2>/dev/null || echo unknown)
    elapsed=""
    start_epoch=\$(date -d '${run_start}' +%s 2>/dev/null || echo "")
    if [ -n "\$start_epoch" ]; then
      secs=\$(( \$(date -u +%s) - start_epoch ))
      elapsed=\$(printf '%dh %02dm' \$((secs/3600)) \$(((secs%3600)/60)))
    fi
    echo " RNAcentral QC report — ${stage}"
    echo " Generated (UTC): \$(date -u '+%Y-%m-%d %H:%M:%S')"
    echo " Loaded into: ${prefix}\$target"
    [ -n "\$elapsed" ] && echo " Elapsed: \$elapsed\""""
}

// Generic single-shot QC report (precompute, genes). errorStrategy 'ignore' so
// QC never breaks the pipeline.
process qc_report {
  tag { stage }
  cache false
  errorStrategy 'ignore'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  publishDir "${params.qc?.publish_dir ?: 'qc-reports'}", mode: 'copy', overwrite: true

  input:
  tuple val(stage), path(query)

  output:
  path("${stage}-qc-report.txt")

  script:
  def run_start = workflow.start ?: ''
  def release = params.get('release', '')
  def notify = params.get('notify', false)
  def fence = '```'
  def slack = notify ? """
  {
    echo "*RNAcentral QC report — ${stage} stage*"
    echo '${fence}'
    cat "\$report"
    echo '${fence}'
  } > slack-report.txt
  rnac notify file slack-report.txt || echo "Slack QC notification failed (non-fatal)"
  """ : "true"
  """
  report="${stage}-qc-report.txt"

  {
    ${qc_header("${stage} stage")}
    psql -q -P pager=off -v run_start='${run_start}' -v release='${release}' -f $query "\$PGDATABASE" 2>&1
  } | tee "\$report"

  ${slack}
  """
}

// Import QC: per-database report (see import.sql). Passes :run_dbs and
// :run_start, so it can't use the generic qc_report process.
process qc_import_report {
  cache false
  errorStrategy 'ignore'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  publishDir "${params.qc?.publish_dir ?: 'qc-reports'}", mode: 'copy', overwrite: true

  input:
  val(_ready)
  path(query)

  output:
  path('import-qc-report.txt')

  script:
  def run_dbs = run_dbs_arg()
  def run_start = workflow.start ?: ''
  def notify = params.get('notify', false)
  def fence = '```'
  def slack = notify ? """
  {
    echo "*RNAcentral QC report — import stage*"
    echo '${fence}'
    cat "\$report"
    echo '${fence}'
  } > slack-report.txt
  rnac notify file slack-report.txt || echo "Slack QC notification failed (non-fatal)"
  """ : "true"
  """
  report="import-qc-report.txt"

  {
    ${qc_header('import stage')}
    echo " Databases imported this run: ${run_dbs ?: '(none detected)'}"
    echo
    psql -q -P pager=off -v run_dbs='${run_dbs}' -v run_start='${run_start}' -f ${query} "\$PGDATABASE" 2>&1
  } | tee "\$report"

  ${slack}
  """
}

// Per-stage entry points.
workflow qc_import {
  take: ready
  main:
    qc_import_report(ready.collect(), file('files/qc/import.sql')) \
    | set { report }
    print_report(report)
  emit: report
}

// Analyze QC is two-phase: snapshot untracked result-table counts at analyze
// start, diff at the end (tracked analyses use pipeline_tracking last_run).
process qc_analyze_before {
  cache false
  errorStrategy 'ignore'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  publishDir "${params.qc?.publish_dir ?: 'qc-reports'}", mode: 'copy', overwrite: true

  input:
  val(_flag)
  path(query)

  output:
  path('analyze-before.csv')

  script:
  """
  psql -q -v snapshot=1 -f ${query} "\$PGDATABASE" 2>&1 \
    | grep -E ',[0-9]' > analyze-before.csv || true
  """
}

process qc_analyze_after {
  cache false
  errorStrategy 'ignore'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  publishDir "${params.qc?.publish_dir ?: 'qc-reports'}", mode: 'copy', overwrite: true

  input:
  val(_ready)
  path('analyze-before.csv')
  path(query)

  output:
  path('analyze-qc-report.txt')

  script:
  def run_start = workflow.start ?: ''
  def notify = params.get('notify', false)
  def fence = '```'
  def slack = notify ? """
  {
    echo "*RNAcentral QC report — analyze stage*"
    echo '${fence}'
    cat "\$report"
    echo '${fence}'
  } > slack-report.txt
  rnac notify file slack-report.txt || echo "Slack QC notification failed (non-fatal)"
  """ : "true"
  """
  report="analyze-qc-report.txt"

  {
    ${qc_header('analyze stage')}
    psql -q -P pager=off -v run_start='${run_start}' -f ${query} "\$PGDATABASE" 2>&1
  } | tee "\$report"

  ${slack}
  """
}

workflow qc_analyze {
  take:
    start_flag   // analyze start  -> before snapshot
    end_flag     // analyses done  -> after snapshot + report
  main:
    qc_analyze_before(start_flag, file('files/qc/analyze.sql')) \
    | set { before }

    qc_analyze_after(end_flag.collect(), before, file('files/qc/analyze.sql')) \
    | set { report }

    print_report(report)
  emit: report
}

workflow qc_precompute {
  take: ready
  main:
    ready \
    | collect \
    | map { ['precompute', file('files/qc/precompute.sql')] } \
    | qc_report \
    | set { report }
    print_report(report)
  emit: report
}

workflow qc_genes {
  take: ready
  main:
    ready \
    | collect \
    | map { ['genes', file('files/qc/genes.sql')] } \
    | qc_report \
    | set { report }
    print_report(report)
  emit: report
}

// Export QC: validate published FTP files -> report (Slack + console); full
// release summary -> file only (too big for Slack).
process qc_export_report {
  cache false
  errorStrategy 'ignore'
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  publishDir "${params.qc?.publish_dir ?: 'qc-reports'}", mode: 'copy', overwrite: true

  input:
  val(_ready)
  path(query)

  output:
  path('export-qc-report.txt'), emit: report
  path('export-summary.txt'),   emit: summary

  script:
  def pubdir = params.export?.ftp?.publish ?: ''
  def run_start = workflow.start ?: ''
  def notify = params.get('notify', false)
  def fence = '```'
  def slack = notify ? """
  {
    echo "*RNAcentral QC report — export stage*"
    echo '${fence}'
    cat "\$report"
    echo '${fence}'
  } > slack-report.txt
  rnac notify file slack-report.txt || echo "Slack QC notification failed (non-fatal)"
  """ : "true"
  """
  report="export-qc-report.txt"
  summary="export-summary.txt"
  pubdir="${pubdir}"

  # 1. Validate the published FTP files -> report (small, Slack + console).
  {
    ${qc_header('export stage')}
    echo " Publish dir: \$pubdir"
    echo
    echo '### Exported files'
    had_problem=0
    for f in \\
      release_notes.txt \\
      sequences/rnacentral_active.fasta.gz \\
      sequences/rnacentral_inactive.fasta.gz \\
      id_mapping/id_mapping.tsv.gz \\
      md5/md5.tsv.gz \\
      rfam/rfam_annotations.tsv.gz \\
      go_annotations/rnacentral_rfam_annotations.tsv.gz \\
      gpi/rnacentral.gpi.gz ; do
      p="\$pubdir/\$f"
      if [ ! -e "\$p" ]; then st=MISSING
      elif [ ! -s "\$p" ]; then st=EMPTY
      elif [ "\${f##*.}" = "gz" ] && ! gzip -t "\$p" 2>/dev/null; then st=CORRUPT
      else continue; fi
      had_problem=1
      printf '  %-8s %s\\n' "\$st" "\$f"
    done
    if [ "\$had_problem" = 0 ]; then
      echo '  All expected files present and valid — ok'
    else
      echo
      echo '  CHECK — export incomplete (problem files listed above)'
    fi
  } | tee "\$report"

  # 2. Full release summary -> file only (too big for Slack).
  {
    ${qc_header('release summary')}
    psql -q -P pager=off -v run_start='${run_start}' -f ${query} "\$PGDATABASE" 2>&1
  } > "\$summary"

  ${slack}
  """
}

workflow qc_export {
  take: ready
  main:
    qc_export_report(ready.collect(), file('files/qc/export.sql')) \
    | set { out }
    print_report(out.report)
  emit: out.report
}
