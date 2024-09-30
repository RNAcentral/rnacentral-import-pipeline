nextflow.enable.dsl=2

process change_status_to_pending {
    input:
    val(_flag)

    output:
    val('done')

    script:
    """
    psql -c "UPDATE litscan_job j SET status='pending' FROM ( SELECT DISTINCT job_id FROM litscan_database WHERE name IN $params.rescan_ids_from ) d WHERE j.job_id = d.job_id;" $PGDB_EMBASSY_USER
    """
}

workflow rescan_ids_by_date {
    take: ready
    emit: done
    main:
      change_status_to_pending(ready) | set{ done }
}

workflow {
    rescan_ids_by_date(Channel.from('ready'))
}
