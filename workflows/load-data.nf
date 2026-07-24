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
  maxForks 2
  cpus 4
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(name), path(ctl), path('raw*.csv')

  output:
  val(name)

  script:
  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  """
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
          echo "Running: \$fn/\$script"
          psql -v ON_ERROR_STOP=1 -f \$script "\$PGDATABASE"
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
    | set { imported_names }

    imported_names
      .flatMap { n -> file("files/import-data/pre-release/*__${n.replace('_', '-')}.sql") }
      .filter { f -> f.exists() }
      .toList()
      .set { pre_scripts }

    imported_names
      .flatMap { n -> file("files/import-data/post-release/*__${n.replace('_', '-')}.sql") }
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
