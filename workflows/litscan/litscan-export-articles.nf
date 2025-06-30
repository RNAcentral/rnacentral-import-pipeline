nextflow.enable.dsl=2

process create_xml_files {
    memory '22GB'
    publishDir "$params.litscan_index", mode: 'copy'

    input:
    val(_flag)

    output:
    val('done')

    script:
    """
    rm -fr "$params.litscan_index"/references_*
    litscan-get-articles.py "$PSYCOPG_CONN" $params.litscan_index
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
