nextflow.enable.dsl=2

process create_xml_files {
    memory '18GB'
    publishDir "$params.litscan_index", mode: 'copy'

    input:
    val(_flag)

    output:
    val('done')

    script:
    """
    rm -fr "$params.litscan_index"/references_*
    litscan-get-articles.py "$PGDB_EMBASSY_USER" $params.litscan_index
    curl -X POST -H 'Content-type: application/json' --data '{"text":"LitScan workflow completed"}' $LITSCAN_SLACK_WEBHOOK
    """
}

workflow export_articles {
    take: ready
    emit: done
    main:
      create_xml_files(ready) | set{ done }
}

workflow {
  export_articles()
}
