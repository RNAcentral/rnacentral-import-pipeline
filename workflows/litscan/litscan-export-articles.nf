nextflow.enable.dsl=2

process create_xml_files {
    memory '22GB'
    publishDir "$params.litscan_index", mode: 'copy', overwrite: true

    input:
    val(_flag)

    output:
    path("references_*.xml.gz")

    script:
    """
    rm -fr "$params.litscan_index"/references_*
    litscan-get-articles.py "$PSYCOPG_CONN" .
    """
}

workflow export_articles {
    take: ready
    emit: reference_files
    main:
      create_xml_files(ready) | set{ reference_files }
}

workflow {
  export_articles()
}
