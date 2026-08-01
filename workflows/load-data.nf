process create_load_tables {
  time '1h'
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  file(create)

  output:
  val('done')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $create "\$PGDATABASE"
  """
}

process merge_and_import {
  tag { name }
  memory 9.GB
  maxForks params.import_data.load_max_forks
  cpus 4
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  // Stage pattern uses the active writer_format so files land as raw*.csv or
  // raw*.parquet and can be globbed by the branch below.
  tuple val(name), path(ctl), path("raw*.${params.writer_format}")

  output:
  tuple val(name), path('rows.count')

  script:
  if (params.writer_format == 'parquet') {
    // TODO(phase-1): add a parquet-equivalent of `rnac validate-pgloader` so
    // this branch verifies rows loaded == rows in the parquet files.
    //
    // No --truncate: create_load_tables drops+recreates all load tables once
    // at the start of load_data, so targets are empty when we append. Keeping
    // --truncate would wipe load_rnacentral_all between the short_sequences
    // and long_sequences runs (both map to it).
    //
    // DuckDB otherwise sizes itself from the machine, not the allocation: it
    // would oversubscribe the cpus and can claim more than the cgroup allows.
    // The memory limit leaves a GB for the Python process.
    """
    load-parquet $name 'raw*.parquet' --count-file rows.count \\
      --threads ${task.cpus} \\
      --memory-limit ${Math.max(1, task.memory.toGiga() - 1)}GB
    """
  } else {
    // pgloader reports no usable count, but a csv that reached here was
    // non-empty by the isEmpty filter below, which is the csv path's gate.
    """
    split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
    echo 1 > rows.count
    """
  }
}

process release {
  time '5d'
  maxForks 1
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory  4.GB

  input:
  path(pre_sql)
  path(post_sql)
  path(limits)

  output:
  val('done')

  when: params.get('should_release', false)

  script:
  def should_release = params.should_release
  def pre = file("work/pre-release")
  def post = file("work/post-release")
  """
  set -o pipefail

  run_sql() {
    local fn="\$1"
    if [[ -s "\$fn" ]]; then
      while IFS='' read -r "script" || [[ -n "\$script" ]]; do
        if [[ ! -z "\$script" ]]; then
          echo "TIMING \$(date +%T) start-sql \$script"
          psql -v ON_ERROR_STOP=1 -f \$script "\$PGDATABASE"
        fi
      done < "\$fn"
    fi
  }

  echo "TIMING \$(date +%T) start release-check"
  ${should_release ? '' : '# ' }rnac --log-level info release check $limits
  echo "TIMING \$(date +%T) start pre-release-sql"
  run_sql "${ Utils.write_ordered(pre, pre_sql.inject([]) { a, fn -> a << fn.getName() }) }"
  echo "TIMING \$(date +%T) start release-run"
  ${should_release ? '' : '# ' }rnac --log-level info release run
  echo "TIMING \$(date +%T) start post-release-sql"
  run_sql "${ Utils.write_ordered(post, post_sql.inject([]) { a, fn -> a << fn.getName() }) }"
  echo "TIMING \$(date +%T) start update-stats"
  ${should_release ? '' : '# ' }rnac --log-level info release update-stats
  echo "TIMING \$(date +%T) done"
  """
}

workflow load_data {
  take: parsed
  main:
    channel.fromPath('files/import-data/limits.json') | set { limits }
    channel.fromPath('files/schema/create_load.sql') | set { schema }

    parsed \
    | filter { f -> !f.isEmpty() } \
    | map { f ->
      def name = f.getBaseName()
      def ctl = file("files/import-data/load/${name.replace('_', '-')}.ctl")
      [[name, ctl], f]
    } \
    | filter { entry ->
      def status = entry[0][1].exists()
      if (!status) {
        log.info "Skipping data ${entry[0][1].getBaseName()}"
      }
      status
    } \
    | groupTuple \
    | map { t -> [t[0][0], t[0][1], t[1]] } \
    | combine(create_load_tables(schema)) \
    | map { n, ctl, fs, _ready -> [n, ctl, fs] } \
    | merge_and_import \
    | filter { _n, rows -> rows.text.trim().toInteger() > 0 } \
    | map { n, _rows -> n } \
    | set { imported_names }

    imported_names
      .flatMap { n -> files("files/import-data/pre-release/*__${n.replace('_', '-')}.sql") }
      .filter { f -> f.exists() }
      .toList()
      .set { pre_scripts }

    imported_names
      .flatMap { n -> files("files/import-data/post-release/*__${n.replace('_', '-')}.sql") }
      .mix(
        channel.fromPath([
          'files/import-data/post-release/000__populate_precompute.sql',
          'files/import-data/post-release/999__cleanup.sql',
        ])
      )
      .filter { f -> f.exists() }
      .toList()
      .set { post_scripts }

    release(pre_scripts, post_scripts, limits) \
    | ifEmpty('no release') \
    | set { released }
  emit: released
}
