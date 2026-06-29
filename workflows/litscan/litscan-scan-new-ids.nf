nextflow.enable.dsl=2

/*
 * Scan newly registered RNA IDs against Europe PMC.
 *
 * Takes the channel of registered-IDs files produced by search_new_ids,
 * fans out one SLURM job per ID, waits for all scans to complete,
 * then emits 'done' for the downstream load_organisms stage.
 *
 * Required environment variables (set in local.config):
 *   PSYCOPG_CONN    - PostgreSQL connection URI
 *
 * The ML classifier model (svc_pipeline.pkl) is bundled in the repo at
 * workflows/litscan/svc_pipeline.pkl and is located automatically.
 */

 process search_job {
  time '4h'
  memory '4GB'
  maxForks 1000

  input:
    path id_chunk

  output:
    path "search_results.csv"

  """
  litscan-search-job.py --job-id $id_chunk --output search_results.csv
  """
 }

process scan_job {
    queue 'datamover'
    memory '8GB'

    input:
    path search_results, stageAs: "input_search_results.dat"

    output:
    path "*.csv"

    script:
    """
    LITSCAN_MODEL="${baseDir}/workflows/litscan/svc_pipeline.pkl" \
        litscan-scan-job.py --search-results input_search_results.dat --xml-directory $params.xml_path
    """
}

process load_job {
  input:
    path csv_files

  output:
    val true

  """
  set -euo pipefail
  psql -v ON_ERROR_STOP=1 "$PSYCOPG_CONN" << EOF
  BEGIN;



  -- Load the articles
  CREATE TEMP TABLE load_litscan_article (LIKE litscan_article INCLUDING DEFAULTS);

  \\copy load_litscan_article (pmcid,title,abstract,author,pmid,doi,year,journal,score,cited_by,retracted,rna_related,probability,type) FROM 'litscan_articles.csv' WITH (FORMAT CSV, HEADER false);

  INSERT INTO litscan_article (pmcid,title,abstract,author,pmid,doi,year,journal,score,cited_by,retracted,rna_related,probability,type)
   SELECT pmcid,title,abstract,author,pmid,doi,year,journal,score,cited_by,retracted,rna_related,probability,type
    FROM load_litscan_article
    ON CONFLICT (pmcid) DO NOTHING;

  -- Load the results table first...
  CREATE TEMP TABLE load_litscan_result (LIKE litscan_result INCLUDING DEFAULTS);

  \\copy load_litscan_result (pmcid, job_id, id_in_title, id_in_abstract, id_in_body) FROM 'litscan_results.csv' WITH (FORMAT CSV, HEADER false);

  INSERT INTO litscan_result  (pmcid, job_id, id_in_title, id_in_abstract, id_in_body)
    SELECT pmcid, job_id, id_in_title, id_in_abstract, id_in_body
    FROM load_litscan_result
    ON CONFLICT (pmcid, job_id) DO NOTHING;

  -- Load abstract sentences
  CREATE TEMP TABLE load_litscan_abstract_sentence (pmcid text, job_id text, sentence text);

  \\copy load_litscan_abstract_sentence (pmcid, job_id, sentence) FROM 'litscan_abstract_sentences.csv' WITH (FORMAT CSV, HEADER false);

  INSERT INTO litscan_abstract_sentence (result_id, sentence)
    SELECT
        r.id,
        s.sentence
    FROM load_litscan_abstract_sentence s
    JOIN litscan_result r
      ON s.pmcid = r.pmcid AND s.job_id = r.job_id;


  -- Load body sentences
  CREATE TEMP TABLE load_litscan_body_sentence (pmcid text, job_id text, sentence text, location text);

  \\copy load_litscan_body_sentence (pmcid, job_id, sentence, location) FROM 'litscan_body_sentences.csv' WITH (FORMAT CSV, HEADER false);

  INSERT INTO litscan_body_sentence (result_id, sentence, location)
    SELECT
        r.id,
        s.sentence,
        s.location
    FROM load_litscan_body_sentence s
    JOIN litscan_result r
      ON s.pmcid = r.pmcid AND s.job_id = r.job_id;

  -- Load hit counts
    CREATE TEMP TABLE load_litscan_hit_counts (job_id text, hit_count int);

    \\copy load_litscan_hit_counts (job_id, hit_count) FROM 'litscan_hit_counts.csv' WITH (FORMAT CSV, HEADER false);

    UPDATE litscan_job
    SET hit_count = COALESCE(litscan_job.hit_count, 0) + staging.hit_count
    FROM load_litscan_hit_counts AS staging
    WHERE litscan_job.job_id = staging.job_id;

  -- Update status
   CREATE TEMP TABLE load_litscan_job_status (job_id text, status text);

  \\copy load_litscan_job_status (job_id, status) from 'litscan_job_status.csv' WITH (FORMAT CSV, HEADER false);

    UPDATE litscan_job
    SET status   = staging.status,
        finished = CASE WHEN staging.status = 'success' THEN now()
                        ELSE litscan_job.finished
                   END
    FROM load_litscan_job_status AS staging
    WHERE litscan_job.job_id = staging.job_id;

  COMMIT;
EOF
  """


}

workflow scan_new_ids {
    take: new_ids_files
    emit: done
    main:
      // Split each registered-IDs file into individual job_id values,
      // then run one SLURM scan_job per ID in parallel.
      new_ids_files
        .splitText(by: 2000, file: true)
        .filter { it.size() > 0 }
        | search_job
        | scan_job
        | flatten
        | collectFile { file ->
            [ "${file.name}", file ]
          }
        | toList
        | load_job

      done = load_job.out
}

workflow {
  // Stand-alone entry point for testing: expects a file path as input param.
  // Example: nextflow run litscan-scan-new-ids.nf --ids_file new_jobs.txt
  scan_new_ids(Channel.fromPath(params.ids_file))
}
