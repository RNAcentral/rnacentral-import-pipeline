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
 *   LITSCAN_MODEL   - path to svc_pipeline.pkl (scikit-learn SVC classifier)
 */

process scan_job {
    input:
    val(job_id)

    output:
    val('done')

    script:
    """
    litscan-scan-job.py --job-id "${job_id}"
    """
}

workflow scan_new_ids {
    take: new_ids_files
    emit: done
    main:
      // Split each registered-IDs file into individual job_id values,
      // then run one SLURM scan_job per ID in parallel.
      new_ids_files
        .splitText { it.trim() }
        .filter { it }
        | scan_job

      // Collect all scan completions; emit 'done' even when no IDs were scanned.
      done = scan_job.out
        .collect()
        .ifEmpty([])
        .map { 'done' }
}

workflow {
  // Stand-alone entry point for testing: expects a file path as input param.
  // Example: nextflow run litscan-scan-new-ids.nf --ids_file new_jobs.txt
  scan_new_ids(Channel.fromPath(params.ids_file))
}
