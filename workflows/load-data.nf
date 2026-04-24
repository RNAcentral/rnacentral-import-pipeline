process create_load_tables {
  time '2d'
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  file(create)

  output:
  val('done')

  """
  psql -v ON_ERROR_STOP=1 -f $create "$PGDATABASE"
  """
}

process merge_and_import {
  tag { name }
  memory 9.GB
  maxForks 2
  cpus 4
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  // Stage pattern uses the active writer_format so files land as raw*.csv or
  // raw*.parquet and can be globbed by the branch below.
  tuple val(name), path(ctl), path("raw*.${params.writer_format}")

  output:
  val(name)

  script:
  if (params.writer_format == 'parquet') {
    // TODO(phase-1): add a parquet-equivalent of `rnac validate-pgloader` so
    // this branch verifies rows loaded == rows in the parquet files.
    //
    // No --truncate: create_load_tables drops+recreates all load tables once
    // at the start of load_data, so targets are empty when we append. Keeping
    // --truncate would wipe load_rnacentral_all between the short_sequences
    // and long_sequences runs (both map to it).
    """
    load-parquet $name 'raw*.parquet'
    """
  } else {
    """
    split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
    """
  }
}

process release {
  time '5d'
  maxForks 1
  when { params.should_release }
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory  4.GB

  input:
  path(pre_sql)
  path(post_sql)
  path(limits)

  output:
  val('done')

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
          echo "Running: \$fn/\$script"
          psql -v ON_ERROR_STOP=1 -f \$script "$PGDATABASE"
        fi
      done < "\$fn"
    fi
  }

  ${should_release ? '' : '# ' }rnac release check $limits
  run_sql "${ Utils.write_ordered(pre, pre_sql.inject([]) { a, fn -> a << fn.getName() }) }"
  ${should_release ? '' : '# ' }rnac release run
  run_sql "${ Utils.write_ordered(post, post_sql.inject([]) { a, fn -> a << fn.getName() }) }"
  ${should_release ? '' : '# ' }rnac release update-stats
  """
}

workflow load_data {
  take: parsed
  emit: released
  main:
    Channel.fromPath('files/import-data/limits.json') | set { limits }
    Channel.fromPath('files/schema/create_load.sql') | set { schema }

    parsed \
    | filter { f -> !f.isEmpty() } \
    | map { f ->
      def name = f.getBaseName()
      def ctl = file("files/import-data/load/${name.replace('_', '-')}.ctl")
      [[name, ctl], f]
    } \
    | filter {
      def status = it[0][1].exists()
      if (!status) {
        log.info "Skipping data ${it[0][1].getBaseName()}"
      }
      status
    } \
    | groupTuple \
    | map { [it[0][0], it[0][1], it[1]] } \
    | combine(create_load_tables(schema)) \
    | map { n, ctl, fs, _ -> [n, ctl, fs] } \
    | merge_and_import \
    | set { imported_names }

    imported_names
      .flatMap { n -> file("files/import-data/pre-release/*__${n.replace('_', '-')}.sql") }
      .filter { f -> f.exists() }
      .toList()
      .set { pre_scripts }

    imported_names
      .flatMap { n -> file("files/import-data/post-release/*__${n.replace('_', '-')}.sql") }
      .mix(
        Channel.fromPath([
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
}
