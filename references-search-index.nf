nextflow.enable.dsl=2

process create_xml_files {
    publishDir "$params.litscan_index", mode: 'copy'

    script:
    """
    rm -fr "$params.litscan_index"/references_*
    references-get-articles.py "$PGDB_EMBASSY_USER" $params.litscan_index
    """
}

workflow {
    create_xml_files()
}
