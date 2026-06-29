#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Rescan selected LitScan jobs against Europe PMC.
 *
 * Selects job_ids from litscan_job by either "finished before a cutoff date"
 * (including finished IS NULL) or membership in a list of source databases
 * (regardless of finished), unions the two sets, and drives them through the
 * standard scan + load + export pipeline.
 *
 * Invocation:
 *   nextflow run litscan-rescan.nf \
 *     --rescan_finished_before 2024-06-11 \
 *     --rescan_databases "ensembl,gtrnadb"
 *
 * At least one of the two params is required.
 *
 * Note on --rescan_databases:
 *   The `litscan_database` table is a snapshot rebuilt by export_metadata at
 *   the end of each full pipeline run (see workflows/litscan/metadata/
 *   load-metadata.ctl, which DROPs and recreates the table from
 *   workflows/litscan/results/*.txt — one file per SQL in
 *   workflows/litscan/queries/, filename basename == database name).
 *   So --rescan_databases can only target names that were present in the
 *   most recent export_metadata run. To add a brand-new source database:
 *     1. Add workflows/litscan/queries/<dbname>.sql.
 *     2. Run the full pipeline (litscan.nf) — this registers and scans any
 *        genuinely new IDs and rebuilds litscan_database including <dbname>.
 *   Rescanning via this workflow does not register new IDs; it only
 *   re-scans existing litscan_job rows.
 *
 * Required environment variables (set in local.config):
 *   PGDATABASE      - psql connection string used by the selector
 *   PSYCOPG_CONN    - PostgreSQL URI used downstream by scan/load
 */

include { scan_new_ids }            from './workflows/litscan/litscan-scan-new-ids'
include { load_organisms }          from './workflows/litscan/litscan-load-organisms'
include { find_organisms }          from './workflows/litscan/litscan-find-organisms'
include { find_retracted_articles } from './workflows/litscan/litscan-retracted-articles'
include { find_manually_annotated } from './workflows/litscan/litscan-manually-annotated'
include { export_articles }         from './workflows/litscan/litscan-export-articles'
include { export_metadata }         from './workflows/litscan/litscan-export-metadata'

params.rescan_finished_before = null
params.rescan_databases       = null

process select_rescan_ids {
  publishDir "$baseDir/workflows/litscan/submit/", mode: 'copy'

  output:
    path 'rescan_ids.txt'

  script:
    def dateClause = params.rescan_finished_before
        ? "(j.finished < '${params.rescan_finished_before}'::timestamp OR j.finished IS NULL)"
        : "FALSE"
    def dbList = params.rescan_databases
        ? params.rescan_databases.split(',').collect { "'${it.trim().replace("'", "''")}'" }.join(',')
        : null
    def dbClause = dbList ? "d.name = ANY(ARRAY[${dbList}])" : "FALSE"
    """
    set -euo pipefail
    psql -v ON_ERROR_STOP=1 -c "COPY (
      SELECT DISTINCT j.job_id
      FROM litscan_job j
      LEFT JOIN litscan_database d ON d.job_id = j.job_id
      WHERE ${dateClause}
         OR ${dbClause}
      ORDER BY 1
    ) TO STDOUT;" \$PGDATABASE > rescan_ids.txt

    count=\$(wc -l < rescan_ids.txt)
    echo "Selected \${count} job_ids for rescan"
    if [ "\${count}" -eq 0 ]; then
      echo "No job_ids matched the rescan criteria; aborting." >&2
      exit 1
    fi
    """
}

workflow {
  if (!params.rescan_finished_before && !params.rescan_databases) {
    error "Must supply at least one of --rescan_finished_before YYYY-MM-DD or --rescan_databases db1,db2"
  }

  select_rescan_ids()
  scan_new_ids(select_rescan_ids.out)
  load_organisms(scan_new_ids.out)
  find_organisms(load_organisms.out)
  find_retracted_articles(find_organisms.out)
  find_manually_annotated(find_retracted_articles.out)
  export_articles(find_manually_annotated.out)
  export_metadata(find_manually_annotated.out, export_articles.out)
}
